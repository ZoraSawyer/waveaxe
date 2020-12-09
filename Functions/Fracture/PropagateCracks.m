function PropagateCracks(dir_cr,nstep,t)
% This function updates the crack geometry, mesh and related data
% structures when the cracks propagate.
%
% INPUTS:
% vel_cr = [ncr x 1] array of crack tip velocities.
% dir_cr = [ncr x 2] array of unit vectors giving the direction of crack
%          propagation
global SMesh % Data for Solid Mesh
global CMesh % Data for crack mesh.

global IOPath % Input/Output files path

%check input data against CMesh data

% if ( size(vel_cr,1) ~= size(dir_cr,1) )
%     disp('Error input array sizes do not match each other ')
%     stop;
% elseif size(vel_cr,1) ~= size(CMesh.tip_nodes,1)
%     disp('Error input array sizes do not match crack mesh data')
%     stop;
% end
% if size(dir_cr,1) ~= size(CMesh.tip_nodes,1)
%     disp('Error input array sizes do not match crack mesh data')
%     stop;
% end

edges = [1,2; 2,3; 3,4; 4,1]; % local Q4 elem connectivity

ntips = size(dir_cr,1);


for i = 1:ntips
    if norm(dir_cr(i,:))
        eID1 = CMesh(i).tip_smesh_e; %enriched element currently containing tip
        tangent = dir_cr(i,:)'/norm(dir_cr(i,:));
        xtip = CMesh(i).nodes(CMesh(i).tip_nodes,:);
        normal = [-tangent(2);tangent(1)];

        % check that propagation direction is +/- 90 degree from previous crack
        % tip direction.
        enodes = SMesh.conn(eID1,:);
        corners = enodes(1:4);
        xI = SMesh.nodes(corners,:);
        gI = SMesh.eLS(eID1,1:4,i);
        if strcmp(SMesh.type,'Q4') || strcmp(SMesh.type,'Q9')
            etp = 'Q4';
        end
        [~,dNdxi] = LagrangeBasis(etp,[0,0],1); % evaluated shape functions at zi
    %     [~,dNdxi]=Q4_LagrangeBasis([0,0]); 

        dxdz = dNdxi'*xI;
        dNdx = dxdz\dNdxi';%(dNdxi*inv(dxdz))'; % Derivative of the shape functions wrt glob coord sys
        gradg = dNdx*gI';

        if gradg'*tangent < 0 
           disp('ERROR in propagatecrack - propagation direction is +/- 90 degree from previous direction');
           stop
        end






        % finding the element edge containing the tip and the neighbour element
        eID2 = 0; % element neighbour to eID1 share edge with the tip.
        gI = SMesh.eLS(eID1,1:4,i);
        fI = SMesh.eLS(eID1,5:8,i);
        enrpos = length(find(gI >= 0));     % number of nodes with positive gI

        for ee = 1:4
            e_sctr = edges(ee,:);
            f1 = fI(e_sctr(1));
            f2 = fI(e_sctr(2));
            g1 = gI(e_sctr(1));
            g2 = gI(e_sctr(2));

            if f1*f2 < 0 
                if enrpos == 3
                    if g1*g2 >= 0
                        eID2 = SMesh.eneighbours(eID1,ee); % neighbouring element ID
                    end
                else
                    if max(g1,g2) >= 0
                        eID2 = SMesh.eneighbours(eID1,ee); % neighbouring element ID
                    end
                end



    %             ee
    %             g1
    %             g2
    %             if ee == 1 || ee == 3
    %                 TEdge_switch = 1;       % a switch indicating whether tip edge is horizontal or 
    %             elseif ee == 2 || ee == 4   % vertical to prevent later propagation along the element edge
    %                 TEdge_switch = 2;
    %             end
            end
        end
        if eID2 <=0
            disp('Error - could not find neighbouring element. ')

            stop;
        end





        % define first approximation of new level sets for eID2
        % in order to find the edge of eID2 which will contain the tip

        nne = 4;
        enodes = SMesh.conn(eID2,:);
        corners = enodes(1:4);
        xtip = xtip + dir_cr(i,:)/norm(dir_cr(i,:));
        for n=1:nne
            xI = SMesh.nodes(corners(n),:);
            r = xI-xtip;
            gI(1,n) = r*tangent;
            fI(1,n) = r*normal;
        end

        % find the edge of element eID2 where the new crack tip will be
        c = 1;
        Xnew = zeros(2,2);
        g = zeros(2,1);
        cut_edges = [0,0];
        enodes = SMesh.conn(eID2,:);
        corners = enodes(1:4);
        xI = SMesh.nodes(corners,:);
        for  ee = 1: 4
            e_sctr = edges(ee,:);
            f1 = fI(e_sctr(1));
            f2 = fI(e_sctr(2));
            g1 = gI(e_sctr(1));
            g2 = gI(e_sctr(2));

            if( f1*f2 < 0 )
                x1 = xI(e_sctr(1),1);
                y1 = xI(e_sctr(1),2);
                x2 = xI(e_sctr(2),1);
                y2 = xI(e_sctr(2),2);
                Mx = (f2-f1)/(x2-x1);
                My = (f2-f1)/(y2-y1);
                Xnew(c,1) = x1-f1/Mx;
                Xnew(c,2) = y1-f1/My;
                g(c) = max(g1,g2);
                cut_edges(c) = ee;
                c = c+1;
            end

        end

        if g(1) > g(2)
            xtip = Xnew(1,:);
            tip_edge = cut_edges(1); % element edge containing the crack tip
        else
            xtip = Xnew(2,:);
            tip_edge = cut_edges(2); % element edge containing the crack tip
        end
        tip_edgenodes = enodes(edges(tip_edge,:));
        tip = CMesh(i).tip_nodes;                                                 % old tip node of the crack
        length_cr = CMesh(i).CrackLength(tip) + norm(xtip-CMesh(i).nodes(tip,:)); % crack length at the new tip node
    %     xtip
        %********************************************
        % element eID2 will now have to be enriched, so set of data structures
        % to make this happend

        % Update list of enriched elements
        % check to make sure that eID2 is not already enriched

        if isempty(find(SMesh.EnrElements == eID2, 1))
            SMesh.EnrElements = [SMesh.EnrElements, eID2];
        else
            disp('Error in propagatecrack - element is already enriched');
            stop;
        end

        % update g level set of eID1 using the new tip location
        enodes = SMesh.conn(eID1,:);
        corners = enodes(1:4);
        for n=1:nne
            xI = SMesh.nodes(corners(n),:);
            r = xI-xtip;
            gI(1,n) = r*tangent;

        end
        SMesh.eLS(eID1,1:4,i) = gI;
        SMesh.nLS(enodes(1:4),1,i) = gI';

        % set the level set of eID2 using the new tip location
        enodes = SMesh.conn(eID2,:);
        corners = enodes(1:4);
        for n=1:nne
            xI = SMesh.nodes(corners(n),:);
            r = xI-xtip;
            gI(1,n) = r*tangent;
            fI(1,n) = r*normal;
        end

        avgsize = (SMesh.eSize(1)+SMesh.eSize(2))/2;    % average element size
        if min(abs(fI)) < .05*avgsize
            [~,index] = min(abs(fI));
            SMesh.EnrExempt(enodes(index)) = 1;
        end

        SMesh.eLS(eID2,1:4,i) = gI;
        SMesh.eLS(eID2,5:8,i) = fI;
        SMesh.nLS(enodes(1:4),1,i) = gI';
        SMesh.nLS(enodes(1:4),2,i) = fI';

        %update the crack mesh by adding a new element and a new node
        old_tip_nID = CMesh(i).tip_nodes;
        new_tip_nID = size(CMesh(i).nodes,1)+1;
        CMesh(i).tip_nodes = new_tip_nID;
        CMesh(i).nodes = [CMesh(i).nodes; xtip];
        CMesh(i).conn = [CMesh(i).conn; [old_tip_nID, new_tip_nID]];
        CMesh(i).GLconn = [CMesh(i).GLconn; new_tip_nID];
        CMesh(i).tip_smesh_e = eID2;
        SMesh.cmesh_e(eID2) = size(CMesh(i).conn,1);
        SMesh.Crnum(eID2) = i;
        CMesh(i).smesh_e = [CMesh(i).smesh_e, eID2];
        CMesh(i).smesh_tipedge = tip_edge;
        CMesh(i).smesh_tipedgenodes = tip_edgenodes;
        CMesh(i).CrackLength = [CMesh(i).CrackLength; length_cr];
        CMesh(i).t0 = [CMesh(i).t0, t];

        CMesh(i).surface_normal = [CMesh(i).surface_normal; normal'];
        SMesh.EnrNodes = [SMesh.EnrNodes, setdiff(enodes(1:4),SMesh.EnrNodes)];
        SMesh.EnrType(enodes(1:4)) = 1;

        CMesh(i).nodes2D = zeros(2*size(CMesh(i).nodes,1),size(CMesh(i).nodes,2));
        CMesh(i).nodes2D(1:2:end-1) = CMesh(i).nodes;
        CMesh(i).nodes2D(2:2:end)   = CMesh(i).nodes;

        cr_conn_pos = 2.*CMesh(i).conn -1;  % positive side connectivity 
        cr_conn_neg = 2.*CMesh(i).conn;     % negative side connectivity

        CMesh(i).conn2D = [cr_conn_pos(:,1) cr_conn_neg cr_conn_pos(:,2)];
    end
end

% Write mesh to files
filename = [IOPath 'mesh.vtk.' num2str(nstep)];
description = 'Solid Mesh Data';
scalardata(1).name = 'ID';
scalardata(1).data = 1:size(SMesh.nodes,1);
scalardata(2).name = 'EnrType';
scalardata(2).data = SMesh.EnrType;
count = 3;
for nc = 1:ntips
    scalardata(count).name = ['gI', num2str(nc)];
    scalardata(count).data = SMesh.nLS(:,1,nc);
    count = count + 1;
    scalardata(count).name = ['fI', num2str(nc)];
    scalardata(count).data = SMesh.nLS(:,2,nc);
    count = count + 1;
end
cell_data(1).name = 'Inclusion';
cell_data(1).data = SMesh.einc;
WriteMesh2VTK(filename,description, SMesh.nodes,SMesh.conn,scalardata,cell_data);

for i = 1:ntips
    filestring = ['crack', num2str(i), '.vtk.', num2str(nstep)];
    filename = [IOPath filestring];
    description = 'Crack Mesh Data';
    S = struct('name',{},'data',{});
    S(1).name = 'ID';
    S(1).data = 1:size(CMesh(i).nodes,1);
    WriteMesh2VTK(filename,description, CMesh(i).nodes,CMesh(i).conn,S);
end

end

