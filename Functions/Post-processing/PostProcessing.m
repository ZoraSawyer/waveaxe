function [SMesh, CMesh, Pvar, Pvar0, Pvar_1, s_dof, f_dof, fdof, enrDOFs, prop] =...
    PostProcessing(Pvar, Pvar0, Pvar_1, stdDOFs, enrDOFs, Q, Q_avg,...
    NR, t, n, save_on, dynamic_ON, dt, SMesh, CMesh, Material, Domain, Control)
% POST PROCESSING Computes stress in solid, fluid flux
% Checks for fracture propagation and updates fracture path if necessary
% Writes data to output files

% input variables:
%
% Pvar: primary variable, Pvar = [d; p1; p2; ... pn] (for n cracks)
% stdDOFs: Standard DOFs
% enrDOFs: Enriched DOFs

% prop: A boolean variable indicating whether the fracture is updated
%       prop = 1, fracture propagation
%       prop = 0, no fracture propagation

% Written by Matin Parchei Esfahani, University of Waterloo, April 2016

disp(['\t', num2str(toc), ': Post processing']);

s_dof  = stdDOFs + enrDOFs;                 % solid DoFs
ncrack = size(CMesh,2);                     % number of cracks
fdof   = zeros(ncrack,1);                   % number of fluid dofs of each crack
f_dof  = 0;                     
for nc = 1:ncrack
    fdof(nc) = size(CMesh(nc).nodes,1);     % number of fluid DoFs for each crack
    f_dof = f_dof + fdof(nc);               % total number of fluid dofs
end
nsd = size(SMesh.nodes,2);                  % number of space dimensions

% ========================= COMPUTE SOLID STRESS ==========================
% Compute Nodal Stresses
disp(['\t', num2str(toc),': Computing stress in solid'])
S = ComputeNodalStress(Pvar(1:s_dof), SMesh, Material, Domain);   % Compute stress at nodes

% in situ stress
Sx  = Domain.InsituStress.Sx;   % in situ stress in x-direction
Sy  = Domain.InsituStress.Sy;   % in situ stress in y-direction
Sxy = Domain.InsituStress.Sxy;  % in situ shear stress (zero on principal planes)

S0 = [Sx; Sy; Sxy];     % in situ stress matrix (voigt)

S_tot = S + repmat(S0,1,size(S,2));   % total stress (S + S0)

% ======================== FRACTURE PROPAGATION ===========================
% Check for fracture propagation
disp(['\t', num2str(toc),': Evaluating fracture propagation criterion'])
if ~isempty(CMesh(1).conn)   % at least one fracture exists
    % CHECKING FOR FRACTURE PROPAGATION
    gplot = 0;  % plot mode is on when gplot = 1; used for debugging purposes
    prop_dir = PropCriterion(S_tot, gplot, SMesh, CMesh, Material);    % Direction of fracture propagation
else                         % no fractrue exists
    prop_dir = [0 0];
end

% Update parameters
if norm(prop_dir) % fracture propagation
    
    % saving old values
    s_dof_old  = s_dof;
    fdof_old   = fdof;
    Pvar_old   = Pvar;
    Pvar0_old  = Pvar0;
    Pvar_1_old = Pvar_1;
    
    [SMesh, CMesh] = PropagateCracks(prop_dir, n, t, SMesh, CMesh, Control.OutPath);       % Updating levelsets and enriched DOFs
                                           % based on the new fracture configuration

    enrDOFs = nsd*length(SMesh.EnrNodes);  % updating number of enriched DOFs
    s_dof   = stdDOFs + enrDOFs;           % updating total number of solid DOFs
    
    % updating number of fluid DOFs
    fdof   = zeros(ncrack,1);              % number of fluid dofs of each crack
    f_dof  = 0;                     
    for nc = 1:ncrack
        fdof(nc) = size(CMesh(nc).nodes,1);     % updating number of fluid DoFs for each crack
        f_dof = f_dof + fdof(nc);               % updating total number of fluid DOFs
    end
    
    % Redefining vectors based on updated DOFs
    Pvar   = zeros(s_dof+f_dof,1);              % new Pvar = [d; p1; p2; ... pn]
    Pvar0  = zeros(s_dof+f_dof,1);              % new Pvar0
    Pvar_1 = zeros(s_dof+f_dof,1);              % new Pvar_1
    
    Pvar(1:s_dof_old)   = Pvar_old(1:s_dof_old);    % d
    Pvar0(1:s_dof_old)  = Pvar0_old(1:s_dof_old);   % d0
    Pvar_1(1:s_dof_old) = Pvar_1_old(1:s_dof_old);  % d_1
    
    count     = 0;
    count_old = 0;
    for nc = 1:ncrack
        Pvar(s_dof+count+1 : s_dof+count+fdof_old(nc)) = ...
        Pvar_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));    % p
        Pvar0(s_dof+count+1 : s_dof+count+fdof_old(nc)) = ...
        Pvar0_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));   % p0
        Pvar_1(s_dof+count+1 : s_dof+count+fdof_old(nc)) = ...
        Pvar_1_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));  % p_1
    
        % padding fluid pressure in the new crack segment with the value of
        % tip pressure to avoid pressure discontinuity
        Pvar(s_dof+count+fdof_old(nc)+1 : s_dof+count+fdof(nc)) = ...
            Pvar_old(s_dof_old+count_old+fdof_old(nc));
        Pvar0(s_dof+count+fdof_old(nc)+1 : s_dof+count+fdof(nc)) = ...
            Pvar0_old(s_dof_old+count_old+fdof_old(nc));
        Pvar_1(s_dof+count+fdof_old(nc)+1 : s_dof+count+fdof(nc)) = ...
            Pvar_1_old(s_dof_old+count_old+fdof_old(nc));
        
        count = count + fdof(nc);
        count_old = count_old + fdof_old(nc);
    end
    prop = 1;   % crack is updated
else
    prop = 0;   % crack is not updated
end

if ~prop && ~save_on
    save_on = 1;
end

% ====================== COMPUTE FRACTURE APERTURE ========================
% Compute fracture aperture
if ~isempty(CMesh(1).conn)     % at least one fracture exists
    disp(['\t', num2str(toc),': Computing fracture aperture'])
    % Computing fracture aperture
    CMesh = Aperture(Pvar(1:s_dof), SMesh, CMesh);
    
    % Find location of the physical tip
    phys_tip = FindPhysicalTip(CMesh, Material);
end

% ================== WRITING RESULTS TO OUTPUT FILES ======================
if save_on && ~mod(n,Control.Postprocessing.OutputFreq)
    
    disp(['\t', num2str(toc),': Saving results to output files'])
    % Writing results of the current time step to file
    dMagCoef = Control.MagCoef;  % displacement magnification coeficient
    xdofs = 1:2:stdDOFs-1;       % x-direction dofs
    ydofs = 2:2:stdDOFs;         % y-direction dofs
    deformedshape = SMesh.nodes + dMagCoef*[Pvar(xdofs) Pvar(ydofs)];      % Compute deformed shape
           

    scalardata(1).name = 'Ux';       scalardata(1).data = Pvar(xdofs);    
    scalardata(2).name = 'Uy';       scalardata(2).data = Pvar(ydofs);
    scalardata(3).name = 'Sxx';      scalardata(3).data = S_tot(1,:);
    scalardata(4).name = 'Syy';      scalardata(4).data = S_tot(2,:);
    scalardata(5).name = 'Sxy';      scalardata(5).data = S_tot(3,:);
    
    % Compute velocity and acceleration vectors (Dynamic analysis)
    if dynamic_ON
        v = 1/dt*(Pvar(1:s_dof) - Pvar0(1:s_dof));
        a = 1/dt^2*(Pvar(1:s_dof) - 2*Pvar0(1:s_dof) + Pvar_1(1:s_dof));
    else
        v = zeros(s_dof,1);
        a = zeros(s_dof,1);
    end
        
    scalardata(6).name = 'Vx';       scalardata(6).data = v(xdofs);    
    scalardata(7).name = 'Vy';       scalardata(7).data = v(ydofs);
    scalardata(8).name = 'ax';       scalardata(8).data = a(xdofs);    
    scalardata(9).name = 'ay';       scalardata(9).data = a(ydofs);
    
    filename = [Control.OutPath 'Solution.vtk.' num2str(n)];    description = 'solution';   
    WriteMesh2VTK(filename, description, deformedshape, SMesh.conn, scalardata);
    
    % Saving NR iterations
    filename = [Control.OutPath 'NR.dat'];
    fileID = fopen(filename,'a');
    fprintf(fileID,'%e\n',NR);
    fclose(fileID);

    if ~isempty(CMesh(1).conn)     % Fracture exists
        count = 0;
        for nc = 1:ncrack
            % Writing fracture results to files
            cr_xdofs = 2*(1:size(CMesh(nc).nodes2D,1))'-1;  % x-direction dofs
            cr_ydofs = 2*(1:size(CMesh(nc).nodes2D,1))';    % y-direction dofs

            dcr = CMesh(nc).dcr;
            w   = CMesh(nc).w;
            deformedshape = CMesh(nc).nodes2D + dMagCoef*[dcr(cr_xdofs) dcr(cr_ydofs)];    % Compute deformed shape

            % nodal values of pressure for the 2D crack mesh
            pressure = zeros(2*fdof(nc),1);

            pressure(1:2:end-1) = Pvar(s_dof+count+1 : s_dof+count+fdof(nc));
            pressure(2:2:end)   = Pvar(s_dof+count+1 : s_dof+count+fdof(nc));

            count = count + fdof(nc);    
            % nodal values of fracture aperture for the 2D crack mesh
            aperture = zeros(2*length(w),1);
            aperture(1:2:end-1) = w;
            aperture(2:2:end)   = w;

            fracdata(1).name = 'Ux';        fracdata(1).data = dcr(cr_xdofs);    
            fracdata(2).name = 'Uy';        fracdata(2).data = dcr(cr_ydofs);
            fracdata(3).name = 'pressure';  fracdata(3).data = pressure;
            fracdata(4).name = 'aperture';  fracdata(4).data = aperture;

            filename = [Control.OutPath 'Fracture' num2str(nc) '.vtk.' num2str(n)];    description = 'fracture';   
            WriteMesh2VTK(filename, description, deformedshape, CMesh(nc).conn2D, fracdata);

            tip_location = CMesh(nc).nodes(CMesh(nc).tip_nodes,:); % Location of the fracture tip

            % Writing location of the fracture tip to file
            tip = CMesh(nc).tip_nodes;                             % tip node of the fracture
            filename = [Control.OutPath 'TipLocation' num2str(nc) '.dat'];
            fileID = fopen(filename,'a');
            fprintf(fileID,'%f %f %f %f %f\n',t,tip_location,CMesh(nc).CrackLength(tip),phys_tip(nc));
            fclose(fileID);
            
            % Saving fracture volume to file
            filename = [Control.OutPath 'Flowrate' num2str(nc) '.dat'];
            fileID = fopen(filename,'a');
            fprintf(fileID,'%f %f %f\n',t,Q(nc),Q_avg(nc));
            fclose(fileID);
            
        end
    end
end
% =========================================================================
end

