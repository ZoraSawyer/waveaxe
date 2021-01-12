function [SMesh, CMesh] = BuildMesh(Mesh, Domain, Control)
% BUILDMESH Builds a structured mesh
% The mesh date is written to a file in VTK form to be read by paraview.
% Last Modified Nov 7, 2012
% copyright Robert Gracie, 2012.
%
% VARIABLE DEFINITIONS:
% nodes_cr = crack mesh nodes
% conn_cr = crack mesh element connectivity
% tip_nodes = [ncr x1 ] list of the indexes into nodes_cr of nodes which 
%             are at crack tips.
% nodes = solid mesh nodes
% conn = solid mesh element connectivity
% eneighbours = [ne x 4] for each element e, row e is a list of the 
%               elements that share an edges with e.  
% TipElements = [ncr x 1] list of solid elements containing the crack tip

disp('** Running Mesh Generator **');
tic;
disp([num2str(toc),': Reading config file...']);

%% MESH GENERATION
  if strcmp(Mesh.Input, 'Gmsh')
      
      [nodes, conn] = LoadMesh(Mesh.FileName, Mesh.nsd, Control.InPath);
      nn  = size(nodes,1);        % number of nodes
      ne  = size(conn,1);         % number of elements
      nne = size(conn,2);         % number of nodes per element
      nex = [];
      ney = [];
      
      nodeconn = zeros(nn, 4);     % nodal connectivity
      elem_inc = zeros(1, ne);     % indicates the material inclusion each element belongs to
      
      % modify mesh to actual wellbore size.
      Lold = Domain.WB(1).radiusmesh; % radius of well in original mesh
      Lnew = Domain.WB(1).radius; % radius of well in new mesh
      xc(1) = Domain.WB(1).center(1); 
      xc(2) =Domain.WB(1).center(2); % coordinates of the center of the well
      dL = Lold-Lnew;
      for I=1:nn
          xI = nodes(I,:);
          dI = xI - xc;
          r = norm(dI);
          if r <= 1.25 
              dr = dL*(r-Lold)-dL;
              c = dI(1)/r; s = dI(2)/r ; %tangent unit vector from node to center
              r = r+dr; % Update new radial distance to center
              nodes(I,1) = xc(1)+r*c; %update nodal coordinates
              nodes(I,2) = xc(2)+r*s;
          end
      end
      
      L   = max(nodes) - min(nodes);           
      Lx  = L(1);                 % length of the model in x-directions
      Ly  = L(2);                 % length of the model in y-directions
      
      if nne == 4                 % element type
          type = 'Q4';            % 4-node quadrilateral element
      elseif nne == 9
          type = 'Q9';            % 9-node quadrilateral element
      end   
  elseif strcmp(Mesh.Input,'Built-in')
      Lx = Domain.Lx;
      Ly = Domain.Ly;
      nex = Mesh.nex;
      ney = Mesh.ney;
      s0x = Mesh.s0x;
      s0y = Mesh.s0y;
      rx = Mesh.rx;
      ry = Mesh.ry;

      [sx, sy, nex, ney] = MeshSegments(Mesh.type, Mesh.elemtype, ...
            Lx, Ly, nex, ney, s0x, s0y, rx, ry);

      switch elemtype
          case 'Q4'
              nne = 4;            % number of nodes per element
              nx  = nex+1;        % number of nodes in x-direction
              ny  = ney+1;        % number of nodes in y-direction
          case 'Q9'
              nne = 9;            % number of nodes per element
              nx  = 2*nex+1;      % number of nodes in x-direction
              ny  = 2*ney+1;      % number of nodes in y-direction
      end

      ne  = nex*ney;              % number of elements
      nn  = nx*ny;                % total number of nodes

      conn = zeros(ne, nne);       % element connectivity
      nodes = zeros(nn, nsd);      % nodal coordinates
      nodeconn = zeros(nn, 4);     % nodal connectivity
      elem_inc = zeros(1, ne);     % indicates the material inclusion each element belongs to

      disp([num2str(toc),': Computing  nodal locations...']);

      c=1;
      for i=1:ny
          for j=1:nx 
              nodes(c,1) = sx(j);
              nodes(c,2) = sy(i);
              c=c+1;
          end    
      end
      if nne == 4
          type = 'Q4';
          disp([num2str(toc),': Defining Q4 element connectivity...']);
          e=1;
          for i=1:ney
              for j=1:nex
                  ii = (i-1)*nx+j;
                  conn(e,1)=ii;
                  conn(e,2)=1+ii;
                  conn(e,3)=1+ii+nx;
                  conn(e,4)=ii+nx;
                  e=e+1;
              end
          end
      elseif nne == 9
          type = 'Q9';
          disp([num2str(toc),': Defining Q9 element connectivity...']);
          e=1;
          for i=1:ney
              for j=1:nex
                  ii = 2*(i-1)*nx+2*(j-1)+1;
                  conn(e,1)=ii;
                  conn(e,2)=2+ii;
                  conn(e,3)=2+ii+2*nx;
                  conn(e,4)=ii+2*nx;
                  conn(e,5)=1+ii;
                  conn(e,6)=2+ii+nx;
                  conn(e,7)=1+ii+2*nx;
                  conn(e,8)=ii+nx;
                  conn(e,9)=ii+1+nx;
                  e=e+1;
              end
          end
      end
  end

%% Define element area
  eArea = zeros(ne,1);    % area of each element

  disp([num2str(toc),': Defining elements of each node...']);
  % list of elements connected to each node
  temp = zeros(nn,1); %temp counter
  for e=1:ne
      enodes = conn(e,:);
      for n=1:nne
          nID = enodes(n);
          c = temp(nID);
          c = c+1;
          temp(nID) = c;
          nodeconn(nID,c) = e;
      end
      eArea(e) = polyarea(nodes(enodes,1),nodes(enodes,2));   % calculate area of each element
  end

%% Define element neighbours
  eneighbours = zeros(ne,4);  % element neighbours (share an edge)
  for e=1:ne
      elist = nodeconn(conn(e,:),:);      
      nnei  = size(elist);
      nnei  = nnei(1)*nnei(2);
      elist = reshape(elist,[nnei,1]);       
      elist = setdiff(unique(elist),[e,0]);  % list of elements which share at least one node with element e 

      pattern = [1,2,3,4,1];
      for i = 1:4
          [r1,~] = find(conn(elist,:) == conn(e,pattern(i)));
          [r2,~] = find(conn(elist,:) == conn(e,pattern(i+1)));

          if isempty(intersect(r1,r2))
              eneighbours(e,i) = 0;
          else
              eneighbours(e,i) = elist(intersect(r1,r2));
          end

      end
  end

%% Define material inclusions
  if Domain.Inclusion_ON
      for n = 1:nn
          for incn = 1:length(Domain.inclusion)   % loop on inclusions
              R = Domain.inclusion(incn).radius;  % radius of the inclusion
              C = Domain.inclusion(incn).center;  % center of the inclusion
              if ((nodes(n,1)-C(1))/R(1))^2 + ((nodes(n,2)-C(2))/R(2))^2 - 1 <= 0 % node located inside an inclusion
                  elem_inc(nodeconn(n,:)) = incn; 
              end
          end
      end
  end

%% Enrichmnet
  % determine nodes to enrich
  if isfield(Domain, 'ncrack')
    ncrack = Domain.ncrack;
  else
    ncrack = 0;
  end
  
  eLS = zeros(ne,2*nne,ncrack);  % level set value for corner nodes of each element [gI,fI] at each node.
  nLS = zeros(nn,2,ncrack);      % level set values at each node for each crack [gI, fI] 
  EnrType = zeros(1,nn);         % Enriched nodes
  EnrExempt = zeros(1,nn);       % Nodes exempted form enrichment
  EnrElements = [];
  EnrNodes    = [];

  Uelements = zeros(ne,1);    % contains crack elements
  eCrnum    = zeros(ne,1);    % indicates the crack number associated to each
                              % element

  if Domain.Fracture_ON  % Model with fracture
          
      disp([num2str(toc),': Computing Level Set values at nodes...']);
      CMesh(ncrack) = struct();
      
      for nc = 1:ncrack       % loop on cracks
          
          TipElements = [];
          nodes_cr    = [];
          temp1       = [];
                  
          normal  = Domain.crack_surface_normal(:,nc);
          tangent = Domain.crack_front_normal(:,nc);
      
          xtip0 = [Domain.x0(nc), Domain.y0(nc)];    % location of the tip
          xtips = [Domain.xs(nc), Domain.ys(nc)];    % location of the start point
          
          fI  = zeros(1,4);   % normal LS
          gI  = zeros(1,4);   % tangent LS
          gI0 = zeros(1,4);   % tangent LS (tip1)
          gIs = zeros(1,4);   % tangent LS (tip2)
          gImax = zeros(1,4); % max tangent LS
          
          switch type
              case 'Q4' 
                  edges = [1,2; 2,3; 3,4; 4,1];
              case 'Q9'
                  edges = [1,2; 2,3; 3,4; 4,1];
          end
         
          count = 0;
          for e=1:ne
             nne = 4;
             enodes = conn(e,:);
             for n=1:nne
                 xI = nodes(enodes(n),:);
                 r1 = xI-xtip0;
                 r2 = xI-xtips;
                 
                 fI(1,n)  = r1*normal;
                 gI0(1,n) = r1*tangent;
                 gIs(1,n) = -r2*tangent;
                 gI(1,n)  = gI0(1,n);
                 gImax(1,n) = max(gI0(1,n),gIs(1,n));
             end
             % set the level set values for this element
             eLS(e,1:2*nne,nc) = [gI,fI];
             nLS(enodes(1:nne),:,nc) = [gI',fI'];

             if max(fI)*min(fI) <0 && min(gImax) < 0
                 
                 if max(gI)*min(gI) < 0 && inpolygon(xtip0(1),xtip0(2),nodes(enodes,1),nodes(enodes,2))
                      TipElements = e;
                 end
                 count = count+1;
                 % enrich this element
                 EnrElements = [EnrElements,e];

                 % define the crack mesh
                 Xnew = zeros(2,2);
                 temp = zeros(2,2);
                 xI = nodes(enodes,:);
                 c = 1;
                 
                 enrpos = length(find(gI >= 0));      % number of nodes with positive gI
                 
                 for ee = 1:4
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
                         temp(c,:) = sort(enodes(e_sctr));   % nodes of the SMesh coresponding to this element edge
                         c = c+1;

                         if ismember(e,TipElements)
                             if enrpos == 3
                                 if g1*g2 >= 0
                                      tip_edge = ee;                        % indicates the edge contains the crack tip
                                      tip_edgenodes = enodes(edges(ee,:));  % indicates nodes of the edge which contains the crack tip
                                 end
                             else
                                  if max(g1,g2) >= 0 % edge containing fracture tip
                                      tip_edge = ee;                        % indicates the edge contains the crack tip
                                      tip_edgenodes = enodes(edges(ee,:));  % indicates nodes of the edge which contains the crack tip
                                  end
                             end
                         end
                     end
                 end
                 if  (Xnew(2,:) - Xnew(1,:))*tangent  > 0
                     nodes_cr = [nodes_cr; Xnew(1,:); Xnew(2,:)];
                     temp1    = [temp1; temp(1,:); temp(2,:)];
                 else
                     nodes_cr = [nodes_cr; Xnew(2,:); Xnew(1,:)];
                     temp1    = [temp1; temp(2,:); temp(1,:)];
                 end
                 
                 if min(abs(fI)) < .05*sqrt(eArea(e)) %.5*(dLx+dLy)
                     [~,index] = min(abs(fI));
                     EnrExempt(enodes(index)) = 1;    % node is exempted from being enriched for stability purposes
                 end
                                
                 EnrNodes = [EnrNodes,setdiff(enodes(1:nne),EnrNodes)];
                 EnrType(enodes(1:nne)) = 1;

                 Uelements(e)   = count;  % starts from 1 for each crack (use eCrnum to find crack number)
                 eCrnum(e)      = nc;     
                 smesh_e(count) = e;
             end

          end

          [temp, index] = unique(temp1,'rows');
          nodes_cr = nodes_cr(index,:);            % crack nodes
          conn_cr = zeros(size(nodes_cr,1)-1,2);   % crack connectivity

          for n = 1:size(nodes_cr,1)   % loop on crack nodes
             for en = 1:size(conn_cr,1)   % loop on crack elements
                 if temp1(2*en-1,:) == temp(n,:)
                     conn_cr(en,1) = n;
                 end
                 if temp1(2*en,:) == temp(n,:)
                     conn_cr(en,2) = n;
                 end
             end
          end

          nodes_cr_2D = zeros(2*size(nodes_cr,1),size(nodes_cr,2));  % 2D crack nodes
          nodes_cr_2D(1:2:end-1,:) = nodes_cr;
          nodes_cr_2D(2:2:end,:)   = nodes_cr;

          cr_conn_pos = 2.*conn_cr -1;  % positive side connectivity 
          cr_conn_neg = 2.*conn_cr;     % negative side connectivity

          conn_cr_2D = [cr_conn_pos(:,1) cr_conn_neg cr_conn_pos(:,2)];  % 2D crack connectivity

          tip_nodes = conn_cr(Uelements(TipElements),2);  % tip nodes

          normalvec = zeros(size(conn_cr,1),2);   % crack normal vectors (negative side)
          for n = 1:size(conn_cr,1)
              normalvec(n,:) = Domain.crack_surface_normal(:,nc);
          end
          
          gI_cr = zeros(size(nodes_cr,1),1);        % tangential LS of the crack nodes
          for nn_cr = 1:size(nodes_cr,1)
              gI_cr(nn_cr) = norm(nodes_cr(nn_cr,:) - nodes_cr(tip_nodes,:));
          end
              
          length_cr = zeros(size(nodes_cr,1),1);    % Fracture length at each node
          [~,glconn] = sort(gI_cr,'descend');
          for cn = 2:size(nodes_cr,1)
              length_cr(glconn(cn)) = length_cr(glconn(cn-1)) + norm(nodes_cr(glconn(cn),:)-nodes_cr(glconn(cn-1),:));
          end
          start_nodes = glconn(1);
                                      
          CMesh(nc).conn2D             = conn_cr_2D;
          CMesh(nc).nodes2D            = nodes_cr_2D;
          CMesh(nc).conn               = conn_cr;
          CMesh(nc).nodes              = nodes_cr;
          CMesh(nc).tip_nodes          = tip_nodes;
          CMesh(nc).start_nodes        = start_nodes;
          CMesh(nc).tip_smesh_e        = TipElements;
          CMesh(nc).smesh_e            = smesh_e;
          CMesh(nc).surface_normal     = normalvec;
          CMesh(nc).smesh_tipedge      = tip_edge;
          CMesh(nc).smesh_tipedgenodes = tip_edgenodes;
          CMesh(nc).CrackLength        = length_cr;
          CMesh(nc).GLconn             = glconn;                      % global connectivity of crack
          CMesh(nc).t0                 = zeros(size(conn_cr,1),1);
          
      end     % loop on cracks
  end

%% Define wellbore nodes
  WBnodes = zeros(nn,1);
  if Domain.Wellbore_ON
     nWB  = length(Domain.WB);
     temp = zeros(nWB,1);
     dWB  = zeros(nn,1);
     for i = 1:nn
         for cWB = 1:nWB
              temp(cWB) = norm(nodes(i,:)-Domain.WB(cWB).center); 
         end
         dWB(i) = min(temp);
     end
     for cWB = 1:nWB
          WBnodes(dWB <= Domain.WB(cWB).radius*1.01) = cWB;
     end
  end

%% Create mesh structure
  SMesh.conn        = conn;
  SMesh.nodes       = nodes;
  SMesh.eneighbours = eneighbours;
  SMesh.nodeconn    = nodeconn;
  SMesh.type        = type;
  SMesh.eLS         = eLS;
  SMesh.nLS         = nLS;
  SMesh.EnrElements = EnrElements;
  SMesh.EnrNodes    = EnrNodes;
  SMesh.EnrType     = EnrType;
  SMesh.EnrExempt   = EnrExempt;
  SMesh.cmesh_e     = Uelements;
  SMesh.eSize       = eArea;%[dLx dLy];
  SMesh.einc        = elem_inc;
  SMesh.Crnum       = eCrnum;
  SMesh.Lx          = Lx;
  SMesh.Ly          = Ly;
  SMesh.nex         = nex;
  SMesh.ney         = ney;
  SMesh.Form        = Mesh.Form;
  SMesh.WBnodes     = WBnodes;

  SMesh = NodeSets(SMesh);

%% Write initial mesh to files
  if ~isfolder(Control.OutPath) 
      mkdir(Control.OutPath)
  end
  filename = [Control.OutPath 'mesh.vtk.0'];
  description = 'Initial Solid Mesh Data';
  scalardata(1).name = 'ID';
  scalardata(1).data = 1:nn;
  scalardata(2).name = 'EnrType';
  scalardata(2).data = EnrType;
  scalardata(3).name = 'Wellbore';
  scalardata(3).data = WBnodes;
  count = 4;
  for nc = 1:ncrack
      scalardata(count).name = ['gI', num2str(nc)];
      scalardata(count).data = nLS(:,1,nc);
      count = count + 1;
      scalardata(count).name = ['fI', num2str(nc)];
      scalardata(count).data = nLS(:,2,nc);
      count = count + 1;
  end
  cell_data(1).name = 'Inclusion';
  cell_data(1).data = elem_inc;
  WriteMesh2VTK(filename, description, SMesh.nodes, SMesh.conn, scalardata, cell_data);

  for nc = 1:ncrack
      filestring = ['crack', num2str(nc), '.vtk.0'];
      filename = [Control.OutPath filestring];
      description = 'Initial Crack Mesh Data';
      S = struct('name',{},'data',{});
      S(1).name = 'ID';
      S(1).data = 1:size(CMesh(nc).nodes,1);
      WriteMesh2VTK(filename,description, CMesh(nc).nodes,CMesh(nc).conn,S);
  end

disp([num2str(toc),': Done.']);