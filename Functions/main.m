clear 
clc
format long
disp('X-FEM HYDRAULIC FRACTURE SIMULATOR ')
tic

%% User input
    ConfigFileName = 'ConfigFile';
 
%% Add directories to search path
    cur_dir = pwd;
    func_dir = cur_dir;
    config_dir = fullfile(cur_dir, '..', 'Config Files');

    addpath(genpath(func_dir));
    addpath(genpath(config_dir));

%% Import config file
    [SMesh, Domain, Material, Control] = feval(ConfigFileName);

%% Build mesh
    disp([num2str(toc),': Loading Mesh and Config file...']);
    [SMesh, CMesh] = BuildMesh(SMesh, Domain, Control);

%% Initialize variables
    % WB normal vectors at each node of the WB
    n_WB = SMesh.nodes(SMesh.WBnodes,:)...
            - repmat(Domain.WB(1).center,length(SMesh.WBnodes),1);
    n_WB = n_WB ./ Domain.WB(1).radius;
    Cwb  = 2*pi*Domain.WB(1).radius / length(SMesh.WBnodes);

    % initializing primary variables
    S0  = [Domain.InsituStress.Sx,  Domain.InsituStress.Sxy;
           Domain.InsituStress.Sxy, Domain.InsituStress.Sy];
    p = zeros(sum([CMesh.fdof]),1);
    count = 0;
    WB_node = zeros(Domain.ncrack,1);
    L       = zeros(Domain.ncrack,1);

    p_hyd = Domain.GravityAcceleration * Domain.Depth * Material.fluid.rho; % Hydro-static pressure

    for n = 1:Domain.ncrack   
        p(count+1:count+CMesh(n).fdof) = p_hyd;     % initial pressure in the n-th fracture
        count = count + CMesh(n).fdof;

        WB_node(n) = CMesh(n).start_nodes;  % Wellbore node of the crack mesh
        tip        = CMesh(n).tip_nodes;    % tip node of the crack mesh
        L(n) = CMesh(n).CrackLength(tip);   % initial length of the fracture
    end

    d = zeros(SMesh.ndof,1);
    p = p_hyd*ones(sum([CMesh.fdof]),1);
    Pvar0  = [d; p];
    Pvar_1 = Pvar0;

    nt     = 0;                               % time step number
    npulse = Control.Wellbore.npulse;         % number of pulses
    iopath = Control.OutPath;

%% Simulation
% Loop through pressure pulses
for np = 1:npulse
    
    fprintf('PULSE %d\n', np);

    % Initialize variables for pulse np
        % change output path to save keep pulse data separate
        new_dir = fullfile(iopath(1:end-1), ['pulse' num2str(np)]);
        if ~isfolder(new_dir)
            mkdir(new_dir)
        end
        Control.OutPath = new_dir;      % Change path

        dt = Control.Time.dtmin;                % time increment
        t0 = -dt;                               % initial time
        t  = t0 + dt;
        
        % set initiation time for leak-off model
        if Domain.Leakoff_ON
            for n = 1:Domain.ncrack
                CMesh(n).t0 = t*ones(1,length(CMesh(n).t0));
            end
        end
        
        save_on = 1;                            % indicates when to save variables
        dynamic_ON = 0;                         % set initial analysis mode to static

        tend = Control.Time.tend;               % final time

        p0t  = zeros(length(t:dt:tend),Domain.ncrack);
        w0t  = zeros(length(t:dt:tend),Domain.ncrack);
        Lt   = zeros(length(t:dt:tend),Domain.ncrack);
        tn   = zeros(length(t:dt:tend),1);

    % Run through timesteps
    while t <= tend

        disp(['\t', num2str(toc),': SOLVING THE SYSTEM OF EQUATIONS'])

        % Initialize variables for timestep nt
            % updating force (NBC)
            Force = zeros(SMesh.ndof + sum([CMesh.fdof]),1);

            % updating fixed DoFs (EBC)
            % tip aperture is zero (restrain enriched dofs)
            % pressure at the wellbore is defined
            tip_count = 0;
            wb_count = 0;
            tipdofs = zeros(1, Domain.ncrack*2*SMesh.nsd);
            pressdofs = zeros(1, Domain.ncrack);

            for n = 1:Domain.ncrack
                tip_sctr = GetScatter(CMesh(n).smesh_tipedgenodes, SMesh);
                % tip dofs of the n-th crack
                tipdofs(tip_count+1:tip_count+2*SMesh.nsd) = tip_sctr(end-3:end);     
                % Wellbore pressure dof of the n-th crack
                pressdofs(n) = SMesh.ndof + wb_count + 1;                  
                tip_count = tip_count + 2*SMesh.nsd;
                wb_count = wb_count + CMesh(n).fdof;
            end    

            fixed_dofs = [SMesh.top_edge_y, SMesh.bot_edge_y, SMesh.left_edge_x, SMesh.right_edge_x, tipdofs, pressdofs];
            
            iter_count  = 0;                        % iteration counter
            converged_Q = 0;                        % flow rate convergence indicator
            max_iter_Q  = 100;                      % maximum number of iterations for flow rate adjustment

        % Update WB pressure
            dpw = Wellbore(t, Material, Control, 'time');      % wellbore pressure at time t
            pw = p_hyd + dpw;                       % update wellbore pressure

            % pressure on well bore
            for i = 1:length(WBnodes)
                fwb = Cwb * (pw*eye(SMesh.nsd) + S0) * n_WB(i,:)';
                Force(WB_x(i)) = fwb(1);
                Force(WB_y(i)) = fwb(2);
            end

            Pvar0(pressdofs) = pw;

        % Newton-Raphson iteration    
            [Pvar, q, NR, NRs, NRf, converged] = NRiter(Force, t, dt, ...
                fixed_dofs, Pvar0, Pvar_1, dynamic_ON, ...
                SMesh, CMesh, Material, Control, Domain);

            if ~converged   % terminates the program if convergence is failed
                return
            end

        % Propagate cracks
            [SMesh, CMesh, Pvar, Pvar0, Pvar_1, prop] = ...
                Propagate(SMesh, CMesh, Material, Domain);

        % Compute fracture aperture
            disp(['\t', num2str(toc),': Computing fracture aperture'])
            CMesh = Aperture(Pvar(1:SMesh.ndof), SMesh, CMesh);

        % Post processing
            % if stopped propagating turn post-processing on
            if ~prop && ~save_on
                save_on = 1;
            end

            if save_on && ~mod(n,Control.Postprocessing.OutputFreq)
                PostProcessing(Pvar, Pvar0, Pvar_1, dynamic_ON, t, n, dt, NR,...
                        SMesh, CMesh, Material, Domain, Control);
            end

        % Move on to the next time step if fracture stopped propagating
            if ~prop

                if save_on
                   for n = 1:Domain.ncrack 
                        p0t(nt+1,n) = Pvar(pressdofs(n));           % wellbore pressure
                        w0t(nt+1,n) = CMesh(n).w(SMesh.WB_node(n));       % wellbore aperture
                        tip         = CMesh(n).tip_nodes;           % tip nod of the crack mesh
                        Lt(nt+1,n)  = CMesh(n).CrackLength(tip);    % fracture length
                        tn(nt+1)    = t;                            % time
                   end
                end

                if Control.Dynamic_ON && ~dynamic_ON               % initial time step of a dynamic analysis
                    Pvar_1 = Pvar;
                    Pvar0  = Pvar;
                else
                    Pvar_1 = Pvar0;
                    Pvar0  = Pvar;
                end

                nt  = nt + 1;
                t0  = t;

                dt = Control.Time.dtmin;
                t = t + dt;
                dynamic_ON = Control.Dynamic_ON;
                save_on = 1;
                fprintf('\nt = %0.6f s \n', t);
            else
                save_on = 0;
            end
    end
end

disp([num2str(toc),': End of analysis'])
