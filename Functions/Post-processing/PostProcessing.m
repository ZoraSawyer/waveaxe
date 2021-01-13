function PostProcessing(Pvar, Pvar0, Pvar_1, dynamic_ON,t, n, dt, NR,...
    SMesh, CMesh, Material, Domain, Control)
% POST PROCESSING Writes data to output files
% input variables:
%
% Pvar: primary variable, Pvar = [d; p1; p2; ... pn] (for n cracks)
% prop: A boolean variable indicating whether the fracture is updated
%       prop = 1, fracture propagation
%       prop = 0, no fracture propagation
% Written by Matin Parchei Esfahani, University of Waterloo, April 2016

    fprintf('\t%.2f: Post processing\n', toc);
   
    % Writing results of the current time step to file
    dMagCoef = Control.MagCoef;  % displacement magnification coeficient
    deformedshape = SMesh.nodes + ...  % Compute deformed shape
                    dMagCoef*[Pvar(SMesh.xdofs) Pvar(SMesh.ydofs)];      
           
    scalardata(1).name = 'Ux';       scalardata(1).data = Pvar(SMesh.xdofs);    
    scalardata(2).name = 'Uy';       scalardata(2).data = Pvar(SMesh.ydofs);
    scalardata(3).name = 'Sxx';      scalardata(3).data = S_tot(1,:);
    scalardata(4).name = 'Syy';      scalardata(4).data = S_tot(2,:);
    scalardata(5).name = 'Sxy';      scalardata(5).data = S_tot(3,:);
    
    % Compute velocity and acceleration vectors (Dynamic analysis)
    if dynamic_ON
        v = 1/dt*(Pvar(1:SMesh.ndof) - Pvar0(1:SMesh.ndof));
        a = 1/dt^2*(Pvar(1:SMesh.ndof) - 2*Pvar0(1:SMesh.ndof) + Pvar_1(1:SMesh.ndof));
    else
        v = zeros(SMesh.ndof,1);
        a = zeros(SMesh.ndof,1);
    end
        
    scalardata(6).name = 'Vx';       scalardata(6).data = v(SMesh.xdofs);    
    scalardata(7).name = 'Vy';       scalardata(7).data = v(SMesh.ydofs);
    scalardata(8).name = 'ax';       scalardata(8).data = a(SMesh.xdofs);    
    scalardata(9).name = 'ay';       scalardata(9).data = a(SMesh.ydofs);
    
    filename = [Control.OutPath 'Solution.vtk.' num2str(n)];    
    description = 'solution';   
    WriteMesh2VTK(filename, description, deformedshape, SMesh.conn, scalardata);
    
    % Saving NR iterations
    filename = [Control.OutPath 'NR.dat'];
    fileID = fopen(filename,'a');
    fprintf(fileID,'%e\n',NR);
    fclose(fileID);

    count = 0;
    for nc = 1:ncrack

        % Writing fracture results to files
        cr_xdofs = 2*(1:size(CMesh(nc).nodes2D,1))'-1;  % x-direction dofs
        cr_ydofs = 2*(1:size(CMesh(nc).nodes2D,1))';    % y-direction dofs

        dcr = CMesh(nc).dcr;
        w   = CMesh(nc).w;
        deformedshape = CMesh(nc).nodes2D + dMagCoef*[dcr(cr_xdofs) dcr(cr_ydofs)];    % Compute deformed shape

        % nodal values of pressure for the 2D crack mesh
        pressure = zeros(2*CMesh(nc).fdof,1);

        pressure(1:2:end-1) = Pvar(SMesh.ndof+count+1 : SMesh.ndof+count+CMesh(nc).fdof);
        pressure(2:2:end)   = Pvar(SMesh.ndof+count+1 : SMesh.ndof+count+CMesh(nc).fdof);

        count = count + CMesh(nc).fdof; 

        % nodal values of fracture aperture for the 2D crack mesh
        aperture = zeros(2*length(w),1);
        aperture(1:2:end-1) = w;
        aperture(2:2:end)   = w;

        fracdata(1).name = 'Ux';        fracdata(1).data = dcr(cr_xdofs);    
        fracdata(2).name = 'Uy';        fracdata(2).data = dcr(cr_ydofs);
        fracdata(3).name = 'pressure';  fracdata(3).data = pressure;
        fracdata(4).name = 'aperture';  fracdata(4).data = aperture;

        filename = [Control.OutPath 'Fracture' num2str(nc) '.vtk.' num2str(n)];   
        description = 'fracture';   
        WriteMesh2VTK(filename, description, deformedshape, CMesh(nc).conn2D, fracdata);

        tip_location = CMesh(nc).nodes(CMesh(nc).tip_nodes,:); % Location of the fracture tip

        % Writing location of the fracture tip to file
        tip = CMesh(nc).tip_nodes;   % tip node of the fracture
        % Find location of the physical tip
        phys_tip = FindPhysicalTip(CMesh, Material);
        filename = [Control.OutPath 'TipLocation' num2str(nc) '.dat'];
        fileID = fopen(filename,'a');
        fprintf(fileID, '%f %f %f %f %f\n', t, tip_location, ...
            CMesh(nc).CrackLength(tip),phys_tip(nc));
        fclose(fileID);                            
    end
end

