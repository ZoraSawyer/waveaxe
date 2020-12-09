function [D] = SolidConstitutive(elem)

% Defines the solid cosntitutive model for a given element 
%
%   Input
%           elem : element
%   Output
%           D : Elasticity matrix

% Written by Matin Parchei Esfahani, University of Waterloo, June 2015

global Material SMesh Domain

inc_sw = 0;
if Domain.Inclusion_ON
    for incn = 1:length(Material.solid.inclusion)   % loop on inclusions
        if SMesh.einc(elem) == incn                 % element belongs to one of the defined inclusions
            E  = Material.solid.constitutive.inclusion(incn).constant(1); % E
            nu = Material.solid.constitutive.inclusion(incn).constant(2); % nu
            inc_sw = 1;
            break
        end
    end
end

if inc_sw == 0  % element does not belong to any of the defined inclusions
    E  = Material.solid.constitutive.constant(1);  % E
    nu = Material.solid.constitutive.constant(2);  % nu
end

switch Material.solid.constitutive.state
    case 'PLANESTRAIN'
        D  = E/((1+nu)*(1-2*nu))*[1-nu  nu   0     ;
                                  nu    1-nu 0     ;
                                  0     0    0.5-nu];
    case 'PLAINSTRESS'
        D  = E/(1-nu^2)*[1  nu 0       ;
                         nu 1  0       ;
                         0  0  (1-nu)/2];
end

end

