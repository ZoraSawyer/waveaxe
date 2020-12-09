clear 
clc
format long
 

cur_dir = pwd;
func_dir = cur_dir;
config_dir = fullfile(cur_dir, '..', 'Config Files');

addpath(genpath(func_dir));
addpath(genpath(config_dir));

disp(' X-FEM HYDRAULIC FRACTURE SIMULATOR ')
tic

disp([num2str(toc),': Loading Mesh and Config file...']);
ConfigFileName = 'ConfigFile';
BuildMesh

nsd = size(SMesh.nodes,2);                          % number of spacedimensions

% Edge nodes of the domain
left_nodes  = find(SMesh.nodes(:,1) == 0);          % left nodes
right_nodes = find(SMesh.nodes(:,1) == SMesh.Lx);   % right nodes
bot_nodes   = find(SMesh.nodes(:,2) == 0);          % bottom nodes
top_nodes   = find(SMesh.nodes(:,2) == SMesh.Ly);   % top nodes

% Left edge DoFs
sctr = GetScatter(left_nodes);
left_edge_x = sctr(1:2:end-1);          % x-DoFs
left_edge_y = sctr(2:2:end);            % y-DoFs

% Right edge DoFs
sctr = GetScatter(right_nodes);
right_edge_x = sctr(1:2:end-1);         % x-DoFs
right_edge_y = sctr(2:2:end);           % y-DoFs

% Top edge DoFs
sctr = GetScatter(top_nodes);
top_edge_x = sctr(1:2:end-1);           % x-DoFs
top_edge_y = sctr(2:2:end);             % y-DoFs

% Bottom edge DoFs
sctr = GetScatter(bot_nodes);
bot_edge_x = sctr(1:2:end-1);           % x-DoFs
bot_edge_y = sctr(2:2:end);             % y-DoFs


% Corner nodes of the domain
LL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == 0);           % lower left corner
LR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == 0);           % lower right corner
UR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == SMesh.Ly);    % upper right corner
UL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == SMesh.Ly);    % upper left corner

% lower left corner DoFs
sctr = GetScatter(LL);
LL_x = sctr(1);                         % x-DoF
LL_y = sctr(2);                         % y-DoF

% lower right corner DoFs
sctr = GetScatter(LR);
LR_x = sctr(1);                         % x-DoF
LR_y = sctr(2);                         % y-DoF

% upper right corner DoFs
sctr = GetScatter(UR);
UR_x = sctr(1);                         % x-DoF
UR_y = sctr(2);                         % y-DoF

% upper left corner DoFs
sctr = GetScatter(UL);
UL_x = sctr(1);                         % x-DoF
UL_y = sctr(2);                         % y-DoF

ncrack = size(CMesh,2);                 % number of cracks

% Wellbore DoFs
WBnodes = find(SMesh.WBnodes == 1);
sctr = GetScatter(WBnodes);
sctr = sctr(1:nsd*length(WBnodes));     % standard DoFs
WB_x = sctr(1:2:end-1);                 % x-DoF
WB_y = sctr(2:2:end);                   % y-DoF

% WB normal vectors at each node of the WB
n_WB = SMesh.nodes(WBnodes,:)-repmat(Domain.WB(1).center,length(WBnodes),1);
n_WB = n_WB ./ Domain.WB(1).radius;


stdDOFs = size(SMesh.nodes,1)*nsd;      % number of standard DoFs
enrDOFs = length(SMesh.EnrNodes)*nsd;   % number of enriched DoFs

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
%     ncr = CMesh(n).surface_normal(1,:);
%     R   = [ ncr(2), ncr(1);
%            -ncr(1), ncr(2)];            % Transformation matrix
%     S_theta = R'*S0*R;                  % Transforming the in situ stress to the fracture coordinates
%     Sn = S_theta(2,2);                  % component of stress normal to the fracture
%     p(count+1:count+fdof(n)) = -Sn;     % initial pressure in the n-th fracture

    p(count+1:count+fdof(n)) = p_hyd;     % initial pressure in the n-th fracture
    count = count + fdof(n);

    WB_node(n) = CMesh(n).start_nodes;  % Wellbore node of the crack mesh
    tip        = CMesh(n).tip_nodes;    % tip nod of the crack mesh
    L(n) = CMesh(n).CrackLength(tip);   % initial length of the fracture

%     Sn_max = min(Sn_max,Sn);            % max normal stress (min function is due to negative sign of in situ stresses) 
end

d = zeros(s_dof,1);
Pvar0  = [d; p];
Pvar_1 = Pvar0;

nt     = 0;                               % time step number
npulse = Control.Wellbore.npulse;         % number of pulses
iopath = IOPath;

for np = 1:npulse
    
    % change output path to save keep pulse data separate
    mkdir(iopath(1:end-1), ['pulse' num2str(np)])   % Creat new folder
    IOPath = [iopath 'pulse' num2str(np) '\'];      % Change path


    % dtmin = Control.Time.dtmin;             % minimum time increment
    % [dpw,~,~] = Wellbore(t0,'time');     % wellbore pressure at time t0
    dt = Control.Time.dtmin;                % time increament
    t0 = -dt;                               % initial time
    t  = t0 + dt;
    % Pvar0(s_dof+WB_node) = ...
    %     Pvar0(s_dof+WB_node)+ dpw;          % impose initial WB pressure
    
    % set initiation time for leak-off model
    if Domain.Leakoff_ON
        for n = 1:ncrack
            CMesh(n).t0 = t*ones(1,length(CMesh(n).t0));
        end
global SMesh CMesh Material Domain Control ConfigFileName OutPath
    end
    
    save_on = 1;                            % indicates when to save variables

    E   = Material.solid.constitutive.constant(1);   % elastic Modulus [Pa]
    nu  = Material.solid.constitutive.constant(2);   % Poisson's ratio of the solid [-]
    % mu  = Material.fluid.mu;                         % fluid viscosity [Pa.s]
    Gc  = Material.solid.FractureEnergy;             % Fracture Energy of the solid [J/m^2]
    Rw  = Domain.WB(1).radius;                       % Radius of the wellbore [m]

    % G   = E/(1+nu)/2;                                % shear modulus [Pa]
    E_p = E/(1-nu^2);                                % plain strain modulus [Pa]

    % lambda = .25*.96*( 2*mu*G^3 / (L^2*(1-nu)^3) )^(1/4) * Q_inj^(-3/4);    % initial value of lambda (dp = lambda * dQ)

    beta = abs(S0(1,1) - S0(2,2)) / sqrt(E_p*Gc/Rw); 
    disp(['Solving the problem for beta = ', num2str(beta)])


    % Q_avg_old  = zeros(ncrack,1);           % average flowrate
    % Vinj_old   = zeros(ncrack,1);           % total injected volume
    % q_tip_old  = 0;                         % initial fluid flux at tip
    % qw         = -Q_inj;                    % injection rate at wellbore
    Vinj = zeros(1,ncrack);
    Vcr  = zeros(1,ncrack);

    tend = Control.Time.tend;               % final time
    dynamic_ON = 0;                         % set initial analysis mode to static

    p0t  = zeros(length(t:dt:tend),ncrack);
    w0t  = zeros(length(t:dt:tend),ncrack);
    Lt   = zeros(length(t:dt:tend),ncrack);
    % q0t  = zeros(length(t:dt:tend),1);
    % qnt  = zeros(length(t:dt:tend),1);
    tn   = zeros(length(t:dt:tend),1);
    % mass = zeros(length(t:dt:tend),1);

    while t <= tend

        % updating force (NBC)
        Force = zeros(s_dof + f_dof,1);

    %     WB_node = CMesh(1).start_nodes;         % Wellbore node of the crack mesh
    %     Force(s_dof+WB_node) = qw;              % injection rate at wellbore

        % updating fixed DoFs (EBC)
        count     = 0;
        count2    = 0;
        tipdofs   = zeros(1,ncrack*2*nsd);
        pressdofs = zeros(1,ncrack);

        for n = 1:ncrack
            sctr = GetScatter(CMesh(n).smesh_tipedgenodes);
            tipdofs(count+1:count+2*nsd) = sctr(end-3:end);     % tip dofs of the n-th crack
            pressdofs(n) = s_dof + count2 + 1;                  % Wellbore pressure dof of the n-th crack
            count = count + 2*nsd;
            count2 = count2 + fdof(n);
        end    

        fixed_dofs = [top_edge_y, bot_edge_y, left_edge_x, right_edge_x, tipdofs, pressdofs];

    %     tol         = 1e-2;                     % tolerance of error for lambda
    %     iter_count  = 0;                        % iteration counter
    %     converged_Q = 0;                        % flow rate convergence indicator
    %     max_iter_Q  = 100;                      % maximum number of iterations for flow rate adjustment
    %     NL          = zeros(max_iter_Q,1);      % stores norm of the error
    %     minerror    = .2;                       % minimum allowable relative change in Qerror for each iteration
    %     C_dn        = .6;                       % correction coefficient to decrease lambda (C_dn < 1)
    %     C_up        = 1.5;                      % correction coefficient to increase lambda (C_up > 1)
    %     Qerror_old  = 1;                        % initial value of error

        disp([num2str(toc),': SOLVING THE SYSTEM OF EQUATIONS'])
    %     while ~converged_Q && iter_count <= max_iter_Q 


            % Update WB pressure
            [dpw,~,~] = Wellbore(t,'time');      % wellbore pressure at time t
    %         pw = -Sn_max + dpw;                     % update wellbore pressure
            pw = p_hyd + dpw;                       % update wellbore pressure
            Cwb  = 2*pi*Domain.WB(1).radius / length(WBnodes);

    %         n_WB = ( SMesh.nodes(WBnodes,:) - repmat(Domain.WB(1).center,length(WBnodes),1) ) ./...
    %             Domain.WB(1).radius;
            % pressure on well bore
            for i = 1:length(WBnodes)
                fwb = Cwb * (pw*eye(nsd) + S0) * n_WB(i,:)';
                Force(WB_x(i)) = fwb(1);
                Force(WB_y(i)) = fwb(2);
            end
            Pvar0(pressdofs) = pw;

           
            % Newton-Raphson iteration    
            [Pvar,q,NR,NRs,NRf,converged] = NRiter(Force, t, dt, s_dof, f_dof, fixed_dofs, Pvar0, Pvar_1, dynamic_ON);
            if ~converged   % terminates the program if convergence is failed
                return
            end

            disp([num2str(toc),': Computing fracture aperture'])
            % Compute fracture aperture
            Aperture(Pvar(1:s_dof));
   %%         
   %% 
    %         CMesh(1).w CMesh(2).w Pvar(s_dof+1:s_dof+f_dof)
    %         disp([num2str(toc),': Computing flow rate'])
            % Compute flow rate
    %         [Q, Q_avg] = ComputeFlowrate(t, dt, Q_avg_old); [Q, Q_avg,
    %         Vcr, Vinj] = ComputeFlowrateVolume(t, dt, Q_avg_old,
    %         Vinj_old, Q_inj);

    %         dQ = Q_inj - Q;         % compute error Qerror =
    %         norm(dQ/Q_inj) tip_node = CMesh.tip_nodes; dm = qw + Q +
    %         q(tip_node);

    %         if Qerror < tol
    %             converged_Q = 1; NL(iter_count+1) = Qerror;
    %             disp([num2str(toc),': CONSERVATION OF FLUID MASS:
    %             CONVERGED'])
    %         else
    %             iter_count     = iter_count + 1; NL(iter_count) = Qerror;
    %                         
    %             disp([num2str(toc),': Adjusting wellbore injection rate:
    %             Attempt ', num2str(iter_count)])
    %             
    %             dQerror = (Qerror_old - Qerror) / Qerror Qerror_old =
    %             Qerror;
    %             
    %             if dQerror < 0
    %                 lambda = lambda * C_dn;     % decrease lambda
    %                 disp([num2str(toc),': Coefficient decreased to
    %                 maintain convergence: Lambda = ', num2str(lambda)])
    %             elseif dQerror < minerror
    %                 lambda = lambda * C_up;     % increase lambda
    %                 disp([num2str(toc),': Coefficient increased to
    %                 increase convergence rate: Lambda = ',
    %                 num2str(lambda)]) Qerror_old = 1;
    %             end

    %             dqw = -lambda * dQ; qw  = qw + dqw;                 %
    %             adjust wellbore injection rate Force(s_dof+WB_node) = qw;
    %             % update wellbore flux q_tip_old = q(tip_node);

    %             dpw = lambda * dQ;      % change in wellbore pressure pw
    %             = pw + dpw;         % updated wellbore pressure
    %             Pvar0(s_dof+1) = pw;

    %         end

    %     end

    %     if ~converged
    %         disp([num2str(toc),': FAILED TO SATISFY CONSERVATION OF FLUID
    %         MASS.']) return
    %     end Q_avg_old = Q_avg; Vinj_old  = Vinj;

        %% Post processing
        [Pvar, Pvar0, Pvar_1, s_dof, f_dof, fdof, enrDOFs, prop] = PostProcessing(Pvar, Pvar0, Pvar_1, ...
            stdDOFs, enrDOFs, Vcr, Vinj, NR, t, nt, save_on, dynamic_ON, dt);

        if ~prop

            if save_on
               for n = 1:ncrack 
                   p0t(nt+1,n) = Pvar(pressdofs(n));           % wellbore pressure
                   w0t(nt+1,n) = CMesh(n).w(WB_node(n));       % wellbore aperture
                   tip         = CMesh(n).tip_nodes;           % tip nod of the crack mesh
                   Lt(nt+1,n)  = CMesh(n).CrackLength(tip);    % fracture length
            %        q0t(nt)  = qw;                        % wellbore injection rate
            %        qnt(nt)  = q(tip_node);               % tip leak-off rate
                   tn(nt+1)    = t;                            % time
            %        mass(nt) = dm;
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
        %     L0  = L;
        %     tip = CMesh(1).tip_nodes;          % tip node of the crack
        %     L   = CMesh(1).CrackLength(tip);   % fracture length
        %     [~,t,~] = Wellbore(L);          % updating wellbore pressure and time
        %     lambda  = lambda * sqrt(L0/L);     % updating lamdba (dp = lambda * dQ)
        %     dt = max(t - t0, dtmin)
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
save([iopath 'CrackMesh.mat'],'CMesh')
disp([num2str(toc),': End of analysis'])
% PlotQ
% PlotFracture