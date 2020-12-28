function [tcoh, dtcohdw] = GetCohesive(w, Material)
% GETCOHESIVE Computes cohesive traction and its derivative wrt fracture
% aperture, given the fraqcture aperture at a point.
%
% Input - w : fracture aperture   
%         tcoh : Cohesive traction
%         dtcohdw : derivative of cohesive traction wrt fracture aperture
%         (dt_coh/dw)

% Written by Matin Parchei Esfahani, University of Waterloo

% Cohesive law constants
switch Material.solid.cohesive.rule
    case 'LINEAR' % linear cohesive law
        fu = Material.solid.TensileStrength;    % ultimate cohesive traction (tensile strength)
        Gc = Material.solid.FractureEnergy;     % fracture energy
        wc = 2*Gc/fu;                           % cohesive zone limit fracture aperture
        ww = 0;                                 % weakening fracture aperture
        k1 = 0;
        k2 = -fu/wc; 
        
    case 'BILINEAR' % Bilinear cohesive law
        fu = Material.solid.TensileStrength;    % ultimate cohesive traction (tensile strength)
        Gc = Material.solid.FractureEnergy;     % fracture energy
        a  = Material.solid.wratio;             % ratio of ww/wc
        wc = 2*Gc/fu;                           % cohesive zone limit fracture openingfu = Material.solid.cohesive.constant(1);  % ultimate cohesive traction 
        ww = a*wc;                              % weakening fracture opening
        k1 = fu/ww;
        k2 = -fu/(wc-ww);
end


% Cohesive traction
if w < ww
    tcoh    = k1*w;              % cohesive traction
    dtcohdw = k1;                % dt_coh/dw 
elseif w < wc
    tcoh    = fu + k2*(w-ww);    % cohesive traction
    dtcohdw = k2;                % dt_coh/dw 
else    % out of cohesive zone
    tcoh    = 0;                 % cohesive traction
    dtcohdw = 0;                 % dt_coh/dw 
end

end

