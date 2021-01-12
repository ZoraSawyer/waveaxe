function [nodes, conn] = LoadMesh(meshfile, nsd, InPath)

% Reads mesh data from a mesh file generated by Gmsh
%
% Input 
%       meshfile : name of the mesh file
%       nsd      : number of space dimensions
%
% Output
%       nodes      : coordinates of mesh nodes
%       conn       : element connectivity

% Written by Matin Parchei Esfahani, University of Waterloo, Sep. 2017

filename = [InPath meshfile '.msh'];
fileID = fopen(filename,'r');

s = textscan(fileID, '%s', 'delimiter', '\n');

n_str = find(strcmp(s{1}, '$Nodes'), 1, 'first');       % start of nodes section
e_str = find(strcmp(s{1}, '$Elements'), 1, 'first');    % start of elements section

nnode = str2double( s{1}(n_str+1) );                    % number of nodes
nelem = str2double( s{1}(e_str+1) );                    % number of elements

nodes = zeros(nnode,nsd);

for i = 1:nnode
    temp = s{1}(n_str+1+i);
    temp = sscanf(temp{1}(1,:), '%f');
    nodes(i,:) = temp(2:nsd+1)';                        % nodal coordinates
end

conn = [];

for i = 1:nelem
    temp = s{1}(e_str+1+i);
    temp = sscanf(temp{1}(1,:), '%d');
    elmtyp = temp(2);                                   % equivalent type number in GMSH
    
    switch elmtyp
        case 1          % L2
            nne = 2;    % number of nodes per element
            edg = 1;
        case 8          % L3
            nne = 3;    % number of nodes per element
            edg = 1;
        case 3          % Q4
            nne = 4;    % number of nodes per element
            edg = 0;
        case 10         % Q9
            nne = 9;    % number of nodes per element
            edg = 0;
        case 15         % single node element
            nne = 1;    % number of nodes per element
            edg = 1;    
    end
    
    if ~edg
        conn = [conn; temp(end-nne+1:end)'];
    end
    
end

end