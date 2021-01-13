function SMesh = NodeSets(SMesh)

SMesh.nsd = size(SMesh.nodes,2);

%% Define node sets and DOF sets
    % Edge nodes of the domain
    SMesh.left_nodes  = find(SMesh.nodes(:,1) == 0);          % left nodes
    SMesh.right_nodes = find(SMesh.nodes(:,1) == SMesh.Lx);   % right nodes
    SMesh.bot_nodes   = find(SMesh.nodes(:,2) == 0);          % bottom nodes
    SMesh.top_nodes   = find(SMesh.nodes(:,2) == SMesh.Ly);   % top nodes

    % Left edge DoFs
    sctr = GetScatter(SMesh.left_nodes, SMesh);
    SMesh.left_edge_x = sctr(1:2:end-1);          % x-DoFs
    SMesh.left_edge_y = sctr(2:2:end);            % y-DoFs

    % Right edge DoFs
    sctr = GetScatter(SMesh.right_nodes, SMesh);
    SMesh.right_edge_x = sctr(1:2:end-1);         % x-DoFs
    SMesh.right_edge_y = sctr(2:2:end);           % y-DoFs

    % Top edge DoFs
    sctr = GetScatter(SMesh.top_nodes, SMesh);
    SMesh.top_edge_x = sctr(1:2:end-1);           % x-DoFs
    SMesh.top_edge_y = sctr(2:2:end);             % y-DoFs

    % Bottom edge DoFs
    sctr = GetScatter(SMesh.bot_nodes, SMesh);
    SMesh.bot_edge_x = sctr(1:2:end-1);           % x-DoFs
    SMesh.bot_edge_y = sctr(2:2:end);             % y-DoFs

    % Corner nodes of the domain
    SMesh.LL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == 0);           % lower left corner
    SMesh.LR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == 0);           % lower right corner
    SMesh.UR = find(SMesh.nodes(:,1) == SMesh.Lx & SMesh.nodes(:,2) == SMesh.Ly);    % upper right corner
    SMesh.UL = find(SMesh.nodes(:,1) == 0        & SMesh.nodes(:,2) == SMesh.Ly);    % upper left corner

    % lower left corner DoFs
    sctr = GetScatter(SMesh.LL, SMesh);
    SMesh.LL_x = sctr(1);                         % x-DoF
    SMesh.LL_y = sctr(2);                         % y-DoF

    % lower right corner DoFs
    sctr = GetScatter(SMesh.LR, SMesh);
    SMesh.LR_x = sctr(1);                         % x-DoF
    SMesh.LR_y = sctr(2);                         % y-DoF

    % upper right corner DoFs
    sctr = GetScatter(SMesh.UR, SMesh);
    SMesh.UR_x = sctr(1);                         % x-DoF
    SMesh.UR_y = sctr(2);                         % y-DoF

    % upper left corner DoFs
    sctr = GetScatter(SMesh.UL, SMesh);
    SMesh.UL_x = sctr(1);                         % x-DoF
    SMesh.UL_y = sctr(2);                         % y-DoF

    % Wellbore DoFs
    sctr = GetScatter(SMesh.WBnodes, SMesh);
    sctr = sctr(1:SMesh.nsd*length(SMesh.WBnodes));     % standard DoFs
    SMesh.WB_x = sctr(1:2:end-1);                 % x-DoF
    SMesh.WB_y = sctr(2:2:end);                   % y-DoF

    SMesh.stdDOFs = size(SMesh.nodes,1)*SMesh.nsd;      % number of standard DoFs
    SMesh.enrDOFs = length(SMesh.EnrNodes)*SMesh.nsd;   % number of enriched DoFs
    SMesh.ndof = SMesh.stdDOFs + SMesh.enrDOFs;

    SMesh.xdofs = 1:2:SMesh.stdDOFs-1;       % x-direction dofs
    SMesh.ydofs = 2:2:SMesh.stdDOFs;         % y-direction dofs

end