function [S,e] = ComputeStressatPoint(X,inputfield,e)

%   Computes stress at a given point, X, in the domain
%
%   Input
%           X : global coordinates of the point at which stress is to be
%               computed
%           inputfield : This input argument can be either the nodal displacement
%                        field, d (number of space dimensions X 1 vector), or the nodal stress 
%                        field (number of nodes X 3 matrix). The code discreminates between these 
%                        two types of input by the field size, and uses the former to compute the
%                        actual stress field and the latter to estimate a smoothed (mapped) stress
%                        field.
%           e : element containing point X (if provided)
%
%   Output
%           S : stress at point X (voigt) S = [Sxx
%                                              Syy
%                                              Sxy]
%           e : element containing point X

% Written by Matin Parchei Esfahani, University of Waterloo, Sep. 2015

global SMesh

nne = size(SMesh.conn,2);                   % number of nodes per element

if nargin == 2                              % if element containing point X is not known
    e = FindElement(X);                     % find element containing point X
end

enodes = SMesh.conn(e,:);                   % element connectivity
sctr   = GetScatter(enodes);                % element DOFs
xI     = SMesh.nodes(enodes,:);             % element nodal coordinates

xi = ParentCoordinates(X,SMesh.type,xI);    % global coordinates transfered to 2D parent coordinates

if size(inputfield,2) == 1                                      % input field is the displacement field, d (nDoF X 1).
    % computing actual stress values
    D = SolidConstitutive(e);                                   % gives stress strain relation based on the constitutive law
    B = Bmatrix(xi, xI, enodes, SMesh.EnrType(enodes), ...
        SMesh.eLS(e,nne+1:end), SMesh.type);                    % B matrix at point X
    S = D*B*inputfield(sctr);                                   % S = [Sxx; Syy; Sxy]
    
else % input field is the stress field, (3 X nnodes for 2D or 6 X nnodes for 3D)
     % computing smoothed (mapped) stress field from nodal values of stress
    N = lagrange_basis(SMesh.type,xi,1);                        % shape functions at point X
    S = inputfield(:,enodes)*N;                                 % S = [Sxx; Syy; Sxy]
end

end
