function [N,J] = Nmatrix(xi,xI,conn,enrnodes,fI,etype,sign)

% Computes FEM/XFEM shape functions and the corresponding Jacobian at a
% given point of an element
%
%   Input
%           xi       : point at which shape functions are computed (in parent coordinates)
%           xI       : nodal coordinates of the element
%           conn     : element connectivity
%           enrnodes : list of enriched nodes of the element
%           fI       : nodal values of normal level set
%           etype    : element type
%           sign     : indicator of the side of the crack (1: positive side; -1: negative side)
%
%
%   Output
%           N : Matrix of shape functions
%           J : Jacobian of the element
%
% N = |N1 0 N2 0 ... Nn 0 N1(H-H1)        0 N2(H-H2)        0 ... Nn(H-Hn)        0|
%     |0 N1 0 N2 ... 0 Nn 0        N1(H-H1) 0        N2(H-H2) ... 0        Nn(H-Hn)|

% Written by Matin Parchei Esfahani, University of Waterloo, Nov 2014
% Last modified Oct. 2017 to include higher order elements. 

global SMesh

if nargin < 7
    sign = 0;
end

nne  = size(conn,2);                        % number of nodes per element
nsd  = size(SMesh.nodes,2);                 % number of space dimensions
enrH = length(find(enrnodes == 1));         % Number of nodes enriched by Heaviside function

if any(enrnodes)    % enriched nodes exist
    count = 0;
    switch sign
        case 0
            [Nv,~] = LagrangeBasis('Q4',xi,1);    % calculating shape functions for interpolating LS (linear)
            f = fI(1:4)*Nv;
            if f>=0
                phi = 0.5;
            elseif f<0
                phi = -0.5;
            end
        case 1      % positive side
            phi = 0.5;
        case -1     % negative side
            phi = -0.5;
    end
end

[Nv,dNdxi] = LagrangeBasis(etype,xi,1);    % calculating shape functions
N = zeros(nsd,nsd*(nne+enrH));
J = dNdxi'*xI;

% forming N matrix
for i = 1:nne
    dofrange = ((i-1)*nsd+1):i*nsd;     % local DOF range (for 2D, 2i-1:2i  ;  for 3D, 3i-2:3i  ;  for 1D, i)
    N(:,dofrange) = Nv(i)*eye(nsd);
    
    if enrnodes(i) == 1                 % Heaviside enrichment
        count = count + 1;
        if fI(i)>0
            phiI = 0.5;
        elseif fI(i)<0
            phiI = -0.5;
        end
        dofrange = nne*nsd + (((count-1)*nsd+1):count*nsd);
        N(:,dofrange) = Nv(i)*(phi-phiI)*eye(nsd);
    end
end
end

