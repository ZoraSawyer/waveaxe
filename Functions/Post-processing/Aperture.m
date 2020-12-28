function CMesh = Aperture(d, SMesh, CMesh)
% APERTURE Computes fracture aperture
%
%   Input
%        d : nodal displacement vector

% Written by Matin Parchei Esfahani, University of Waterloo, Sep. 2016

etype = SMesh.type;
ncrack = size(CMesh,2);     % number of cracks

% Compute fracture deformation

% dcr = |dcr_x_pos(1)|  
%       |dcr_y_pos(1)|   
%       |dcr_x_neg(1)|   
%       |dcr_y_neg(1)|    
%       |    ...     |         
%       |dcr_x_pos(n)|       
%       |dcr_y_pos(n)|        
%       |dcr_x_neg(n)|       
%       |dcr_y_neg(n)|       

if strcmp(etype,'Q4')
    fLSrange = 5:8;                             % range of normal LSs for each elements
elseif strcmp(etype,'Q9')
    fLSrange = 5:8;                             % range of normal LSs for each elements
end

for nc = 1:ncrack
    % Memory pre-allocation
    dcr = zeros(2*size(CMesh(nc).nodes2D,1),1); % crack deformation
    w   = zeros(size(CMesh(nc).nodes,1),1);     % crack aperture

    count   = 1;
    count_w = 1;

    for ncr = 1:length(CMesh(nc).smesh_e)
        crnodes = CMesh(nc).conn(ncr,:);        % crack element connectivity
        crnodes = CMesh(nc).nodes(crnodes,:);   % global coordinates of the crack nodes
        emesh   = CMesh(nc).smesh_e(ncr);       % coresponding mesh element
        enodes  = SMesh.conn(emesh,:);          % element nodes
        xI      = SMesh.nodes(enodes,:);        % nodal coordinates of the element (global)
        sctr    = GetScatter(enodes, SMesh);           % element DOFs

        if ncr == length(CMesh(nc).smesh_e)     % element which contains fracture tip
            i = 2;
        else
            i = 1;
        end

        for npt = 1:i
            [xi] = ParentCoordinates(crnodes(npt,:),etype,xI);    % Parent coordinates of the crack node

            [N,~] = Nmatrix(xi,xI,enodes,SMesh.EnrType(enodes),...
                SMesh.eLS(emesh,fLSrange,nc),etype,1);     % positive side
            dcr_pos = N*d(sctr);                           % crack nodes displacement (positive side)
            dcr([2*count-1,2*count]) = dcr_pos;         
            count = count + 1;

            [N,~] = Nmatrix(xi,xI,enodes,SMesh.EnrType(enodes),...
                SMesh.eLS(emesh,fLSrange,nc),etype,-1);    % negative side
            dcr_neg = N*d(sctr);                           % crack nodes displacement (negative side)
            dcr([2*count-1,2*count]) = dcr_neg;
            count = count + 1;

            w(count_w) = CMesh(nc).surface_normal(ncr,:)*(dcr_pos - dcr_neg);   % fracture aperture
            count_w = count_w + 1;
        end
    end
    CMesh(nc).w = w;
    CMesh(nc).dcr = dcr;
end

end

