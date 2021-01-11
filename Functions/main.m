clear 
clc
format long
disp(' X-FEM HYDRAULIC FRACTURE SIMULATOR ')
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

    nsd = size(SMesh.nodes,2);              % number of space dimensions
    ncrack = size(CMesh,2);                 % number of cracks

%% Define node sets and DOF sets
    % Edge nodes of the domain
    left_nodes  = find(SMesh.nodes(:,1) == 0);          % left nodes
    right_nodes = find(SMesh.nodes(:,1) == SMesh.Lx);   % right nodes
    bot_nodes   = find(SMesh.nodes(:,2) == 0);          % bottom nodes
    top_nodes   = find(SMesh.nodes(:,2) == SMesh.Ly);   % top nodes

    % Left edge DoFs
    sctr = GetScatter(left_nodes, SMesh);
    left_edge_x = sctr(1:2:end-1);          % x-DoFs
    left_edge_y = sctr(2:2:end);            % y-DoFs

    % Right edge DoFs
    sctr = GetScatter(right_nodes, SMesh);
    right_edge_x = sctr(1:2:end-1);         % x-DoFs
    right_edge_y = sctr(2:2:end);           % y-DoFs

    % Top edge DoFs
    sctr = GetScatter(top_nodes, SMesh);
    top_edge_x = sctr(1:2:end-1);           % x-DoFs
    top_edge_y = sctr(2:2:end);             % y-DoFs

    % Bottom edge DoFs
    sctr = GetScatter(bot_nodes, SMesh);
    bot_edge_x = sctr(1:2:end-1);           % x-DoFs
    bot_edge_y = sctr(2:2:end);             % y-DoFs


    % Corner nodes of the domain
    LL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == 0);           % lower left corner
    LR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == 0);           % lower right corner
    UR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == SMesh.Ly);    % upper right corner
    UL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == SMesh.Ly);    % upper left corner

    % lower left corner DoFs
    sctr = GetScatter(LL, SMesh);
    LL_x = sctr(1);                         % x-DoF
    LL_y = sctr(2);                         % y-DoF

    % lower right corner DoFs
    sctr = GetScatter(LR, SMesh);
    LR_x = sctr(1);                         % x-DoF
    LR_y = sctr(2);                         % y-DoF

    % upper right corner DoFs
    sctr = GetScatter(UR, SMesh);
    UR_x = sctr(1);                         % x-DoF
    UR_y = sctr(2);                         % y-DoF

    % upper left corner DoFs
    sctr = GetScatter(UL, SMesh);
    UL_x = sctr(1);                         % x-DoF
    UL_y = sctr(2);                         % y-DoF

    % Wellbore DoFs
    WBnodes = find(SMesh.WBnodes == 1);
    sctr = GetScatter(WBnodes, SMesh);
    sctr = sctr(1:nsd*length(WBnodes));     % standard DoFs
    WB_x = sctr(1:2:end-1);                 % x-DoF
    WB_y = sctr(2:2:end);                   % y-DoF

    stdDOFs = size(SMesh.nodes,1)*nsd;      % number of standard DoFs
    enrDOFs = length(SMesh.EnrNodes)*nsd;   % number of enriched DoFs

%% Initialize variables
    % WB normal vectors at each node of the WB
    n_WB = SMesh.nodes(WBnodes,:)-repmat(Domain.WB(1).center,length(WBnodes),1);
    n_WB = n_WB ./ Domain.WB(1).radius;

    s_dof = stdDOFs + enrDOFs;              % total number of solid DoFs
    fdof  = zeros(ncrack,1);                
    f_dof = 0;
    for n = 1:ncrack
        fdof(n) = size(CMesh(n).nodes,1);   % number of dofs of the n-th each crack
        f_dof   = f_dof + fdof(n);          % total number of fluid DoFs
    end

    % initializing primary variables
    S0  = [Domain.InsituStress.Sx,  Domain.InsituStress.Sxy;
           Domain.InsituStress.Sxy, Domain.InsituStress.Sy];
    p = zeros(f_dof,1);
    count = 0;
    WB_node = zeros(ncrack,1);
    L       = zeros(ncrack,1);
    Sn_max  = 0;

    p_hyd = Domain.GravityAcceleration * Domain.Depth * Material.fluid.rho; % Hydro-static pressure

    for n = 1:ncrack   
        p(count+1:count+fdof(n)) = p_hyd;     % initial pressure in the n-th fracture
        count = count + fdof(n);

        WB_node(n) = CMesh(n).start_nodes;  % Wellbore node of the crack mesh
        tip        = CMesh(n).tip_nodes;    % tip nod of the crack mesh
        L(n) = CMesh(n).CrackLength(tip);   % initial length of the fracture
    end

    d = zeros(s_dof,1);
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
        Control.OutPath = [iopath 'pulse' num2str(np) '\'];      % Change path

        dt = Control.Time.dtmin;                % time increament
        t0 = -dt;                               % initial time
        t  = t0 + dt;
        
        % set initiation time for leak-off model
        if Domain.Leakoff_ON
            for n = 1:ncrack
                CMesh(n).t0 = t*ones(1,length(CMesh(n).t0));
            end
        end
        
        save_on = 1;                            % indicates when to save variables

        E   = Material.solid.constitutive.constant(1);   % elastic Modulus [Pa]
        nu  = Material.solid.constitutive.constant(2);   % Poisson's ratio of the solid [-]
        Gc  = Material.solid.FractureEnergy;             % Fracture Energy of the solid [J/m^2]
        Rw  = Domain.WB(1).radius;                       % Radius of the wellbore [m]
        E_p = E/(1-nu^2);                                % plain strain modulus [Pa]

        Vinj = zeros(1,ncrack);
        Vcr  = zeros(1,ncrack);

        tend = Control.Time.tend;               % final time
        dynamic_ON = 0;                         % set initial analysis mode to static

        p0t  = zeros(length(t:dt:tend),ncrack);
        w0t  = zeros(length(t:dt:tend),ncrack);
        Lt   = zeros(length(t:dt:tend),ncrack);
        tn   = zeros(length(t:dt:tend),1);

    % Run through timesteps
    while t <= tend

        disp([num2str(toc),': SOLVING THE SYSTEM OF EQUATIONS'])
        % Initialize variables for timestep nt
            % updating force (NBC)
            Force = zeros(s_dof + f_dof,1);

            % updating fixed DoFs (EBC)
            count     = 0;
            count2    = 0;
            tipdofs   = zeros(1,ncrack*2*nsd);
            pressdofs = zeros(1,ncrack);

            for n = 1:ncrack
                sctr = GetScatter(CMesh(n).smesh_tipedgenodes, SMesh);
                tipdofs(count+1:count+2*nsd) = sctr(end-3:end);     % tip dofs of the n-th crack
                pressdofs(n) = s_dof + count2 + 1;                  % Wellbore pressure dof of the n-th crack
                count = count + 2*nsd;
                count2 = count2 + fdof(n);
            end    

            fixed_dofs = [top_edge_y, bot_edge_y, left_edge_x, right_edge_x, tipdofs, pressdofs];
            
            iter_count  = 0;                        % iteration counter
            converged_Q = 0;                        % flow rate convergence indicator
            max_iter_Q  = 100;                      % maximum number of iterations for flow rate adjustment


        % Update WB pressure
            dpw = Wellbore(t, Material, Control, 'time');      % wellbore pressure at time t

            pw = p_hyd + dpw;                       % update wellbore pressure

            Cwb  = 2*pi*Domain.WB(1).radius / length(WBnodes);

            % pressure on well bore
            for i = 1:length(WBnodes)
                fwb = Cwb * (pw*eye(nsd) + S0) * n_WB(i,:)';
                Force(WB_x(i)) = fwb(1);
                Force(WB_y(i)) = fwb(2);
            end

            Pvar0(pressdofs) = pw;

        % Newton-Raphson iteration    
            [Pvar, q, NR, NRs, NRf, converged] = NRiter(Force, t, dt, ...
                s_dof, f_dof, fixed_dofs, Pvar0, Pvar_1, dynamic_ON, ...
                SMesh, CMesh, Material, Control, Domain);
            if ~converged   % terminates the program if convergence is failed
                return
            end

        disp([num2str(toc),': Computing fracture aperture'])

        % Compute fracture aperture
        CMesh = Aperture(Pvar(1:s_dof), SMesh, CMesh);

        %% Post processing
        [SMesh, CMesh, Pvar, Pvar0, Pvar_1, s_dof, f_dof, fdof, enrDOFs, prop] = PostProcessing(Pvar, Pvar0, Pvar_1, ...
            stdDOFs, enrDOFs, Vcr, Vinj, NR, t, nt, save_on, dynamic_ON, dt, SMesh, CMesh, Material, Domain, Control);

        % move on to the next time step if fracture stopped propagating
        if ~prop

            if save_on
               for n = 1:ncrack 
                    p0t(nt+1,n) = Pvar(pressdofs(n));           % wellbore pressure
                    w0t(nt+1,n) = CMesh(n).w(WB_node(n));       % wellbore aperture
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
