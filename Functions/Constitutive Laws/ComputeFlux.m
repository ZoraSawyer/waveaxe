function [] = ComputeFlux(p)
% Computes fluid flux in a 1D fracture
%
% Input
%       p : fluid pressure vector

% Written by Matin Parchei Esfahani, University of Waterloo, May 2017

global CMesh Material

ncrack = size(CMesh,2);                   % number of cracks
mu   = Material.fluid.mu;                 % viscosity of the fluid
Gc   = Material.solid.FractureEnergy;     % fracture energy
fu   = Material.solid.TensileStrength;    % ultimate cohesive traction (tensile strength)
wmin = 2*Gc/fu;                           % minumum fracture aperture for fluid model (wmin = wc of the cohesive model)
wmin = max(Material.fluid.constitutive.wmin,wmin);

for nc = 1:ncrack
    ne_c = size(CMesh(nc).conn,1);        % number of crack elements
    q = zeros(ne_c,1);                    % flux in each element
    w = max(CMesh(nc).w,wmin);            % fracture aperture
    glc = CMesh.GLconn;                   % global connectivity of the fracture
    
    for e = 1:ne_c
        le = CMesh(nc).CrackLength(glc(e+1)) - CMesh(nc).CrackLength(glc(e));   % length of the crack element
        dpdx = (p(glc(e+1)) - p(glc(e)))/le;                                    % pressure gradient
        Dw   = (w(glc(e))+w(glc(e+1)))^3/(8*12*mu);                             % Dw = w^3/12*mu (using average aperture)
        q(e) = -Dw * dpdx;                                                      % fluid flux in the element
    end
    CMesh(nc).q = q;
end

end

