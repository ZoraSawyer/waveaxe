function [ S ] = ComputeNodalStress( d )

% Returns a smooth matrix of nodal stresses
%
%   Input
%           d : solid displacement vector
%
%   Output
%           S : values of stress at each node of the mesh
%           
%           |Sxx_1 Sxx_2 ... Sxx_n|
%       S = |Syy_1 Syy_2 ... Syy_n|
%           |Sxy_1 Sxy_2 ... Sxy_n|

% Copyright Matin Parchei Esfahani, University of Waterloo, June 2015

global SMesh

nn  = size(SMesh.nodes,1);          % number of nodes
ne  = size(SMesh.conn,1);           % number of elements
nne = size(SMesh.conn,2);           % number of nodes per element

S = zeros(3,nn);
count = zeros(1,nn);

for e = 1:ne                        % loop on elements
    
    enodes = SMesh.conn(e,:);       % element connectivity
    sctr   = GetScatter(enodes);    % element DOFs
    xI     = SMesh.nodes(enodes,:); % element nodal coordinates
    
    if SMesh.Crnum(e)               % if element contains the crack 
        crnum = SMesh.Crnum(e);     % crack number
    else                            % if element doesn't contain the crack
        crnum = 1;                  % LS of the first crack is passed to Bmatrix
    end
    
    D = SolidConstitutive(e);       % gives stress strain relation based on the constitutive law
    
    switch nne
        case 4                      % Q4 element
            xi = [-1 -1;
                   1 -1;
                   1  1;
                  -1  1];
            etype = 'Q4';
            fLSrange = 5:8;         % range of normal LSs for the element
            
        case 9                      % Q9 element
            xi = [-1 -1;
                   1 -1;
                   1  1;
                  -1  1;
                   0 -1;
                   1  0;
                   0  1;
                  -1  0;
                   0  0];
            etype = 'Q9';
            fLSrange = 5:8;         % range of normal LSs for the element
    end
        
    for n = 1:nne
        [B,~] = Bmatrix(xi(n,:), xI, enodes, SMesh.EnrType(enodes),...
            SMesh.eLS(e,fLSrange,crnum), etype);          % B matrix at nodes
        S(:,enodes(n)) = S(:,enodes(n)) + D*(B*d(sctr));  % stress tensor at each node (Voigt)
        count(enodes(n)) = count(enodes(n)) + 1;          % number of elements containing this node
    end
            
end

% Averaging stress over nodes
S(1,:) = S(1,:)./count;
S(2,:) = S(2,:)./count;
S(3,:) = S(3,:)./count;

end

