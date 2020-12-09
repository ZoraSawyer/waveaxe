function [ phys_tip ] = FindPhysicalTip()

% Finds location of the physical tip (t_coh = 0) 
% 
% output : Phys_tip: location of the physical tip in crack coordinate
% system (s=0 at the WB and s=L at the mathematical tip)

% Written by Matin Parchei Esfahani, University of Waterloo, Oct. 2017

global CMesh Material

ncrack   = size(CMesh,2);                     % number of cracks
Gc       = Material.solid.FractureEnergy;     % fracture energy
fu       = Material.solid.TensileStrength;    % ultimate cohesive traction (tensile strength)
wcr      = 2*Gc/fu;                           % critical aperture
phys_tip = zeros(ncrack,1);

for nc = 1:ncrack
    
    sw = 0;
    for i = length(CMesh(nc).w):-1:1
        if CMesh(nc).w(CMesh(nc).GLconn(i)) >= wcr
            sw = 1;
            break
        end
    end 
    
    if i == length(CMesh(nc).w)                 % no cohesive zone
        tip = CMesh(nc).tip_nodes;              % tip node of the crack
        phys_tip(nc) = CMesh(nc).CrackLength(tip);
    elseif ~sw                                  % the whole crack length is cohesive zone
        str = CMesh(nc).start_nodes;            % first node of the crack
        phys_tip(nc) = CMesh(nc).CrackLength(str);
    else                                        
        w1 = CMesh(nc).w(CMesh(nc).GLconn(i));                    
        L1 = CMesh(nc).CrackLength(CMesh(nc).GLconn(i));
        w2 = CMesh(nc).w(CMesh(nc).GLconn(i+1));
        L2 = CMesh(nc).CrackLength(CMesh(nc).GLconn(i+1));
        phys_tip(nc) = L1 + (wcr - w1) * (L2 - L1) / (w2 - w1);     % location of the physical tip (fracture coordinate - linear interpolation)
    end
end
        
end

