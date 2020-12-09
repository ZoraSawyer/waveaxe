function WriteMesh2VTK_v2(filename,description, nodes,conn,scalar_data,cell_data)
% This funciton write FEM mesh data to a file in VTK format
% INPUTS:
% nodes = array of nodal coordinates (n x nsd)
% conn = connectivity matrix (ne x nne)
%
% Last Modified Nov 7, 2012
% copyright Robert Gracie, 2012.

% modified to accommodate higher order elements, Oct 2017, Matin Parchei Esfahani

if nargin < 6
    cell_data = [];
end

%msg = sprintf('Writting results to "%s"\n', filename);
%disp(msg)
% disp([': Writting results to ',filename]);
nn = size(nodes,1); % number of nodes (points)
nsd = size(nodes,2); % number of space dimensions
ne = size(conn,1); % number of elements (cells)
nne = size(conn,2); % number of nodes per element 
nd = length(scalar_data); % number of scalar datum defined at each node (point);
nc = length(cell_data);   % number of scalar datum defined at each cell (element)

fid = fopen(filename, 'w');


fprintf(fid,'%s\n','# vtk DataFile Version 2.0');
fprintf(fid,'%s\n',['HFX mesh description: ',description]);
fprintf(fid,'%s\n','ASCII');

fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');

text = ['POINTS ',num2str(nn),' float'];
fprintf(fid,'%s\n',text);

if nsd == 2
    nodes = [nodes,zeros(nn,1)];
end


fprintf(fid,'%f %f %f\n',nodes');


text = ['CELLS ',num2str(ne), ' ',num2str( (nne+1)*ne)];
fprintf(fid,'%s\n',text);





if nne == 4     % Q4 
  outputformat = '%d %d %d %d %d \n';
%   conn = [4*ones(ne,1),conn(:,1:4)-ones(ne,4)];
  conn = [4*ones(ne,1),conn-ones(ne,4)];
  cell_type = 9;
elseif nne == 9 % Q9
  outputformat = '%d %d %d %d %d %d %d %d %d %d \n';  
  conn = [9*ones(ne,1),conn-ones(ne,9) ];
  cell_type = 23;    
elseif nne == 2 % L2
  outputformat = '%d %d %d  \n';
  conn = [nne*ones(ne,1),conn-ones(ne,nne)];
  cell_type = 3;
else            % B8
  conn = [nne*ones(ne,1),conn-ones(ne,nne)];
  outputformat = '%d %d %d %d %d %d %d %d %d \n';
  cell_type = 12;
end
   
fprintf(fid,outputformat,conn');

text = ['CELL_TYPES ',num2str(ne)];
fprintf(fid,'%s\n',text);

fprintf(fid,'%d\n',cell_type*ones(1,ne));

if nd > 0
text = ['POINT_DATA ',num2str(nn)];
fprintf(fid,'%s\n',text);
end

for i = 1:nd
    text = ['SCALARS ',scalar_data(i).name,' float 1'];
    fprintf(fid,'%s\n',text);
    text = 'LOOKUP_TABLE default';
    fprintf(fid,'%s\n',text);
    fprintf(fid,'%f\n',scalar_data(i).data');
end

if nc > 0
text = ['CELL_DATA ',num2str(ne)];
fprintf(fid,'%s\n',text);
end

for i = 1:nc
    text = ['SCALARS ',cell_data(i).name,' float 1'];
    fprintf(fid,'%s\n',text);
    text = 'LOOKUP_TABLE default';
    fprintf(fid,'%s\n',text);
    fprintf(fid,'%f\n',cell_data(i).data');
end


fclose(fid);

% # vtk DataFile Version 2.0
%     Cube example
%     ASCII