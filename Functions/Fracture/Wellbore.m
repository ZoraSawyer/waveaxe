function [pw, t, Q] = Wellbore(var, Material, Control, mode)
% WELLBORE Returns wellbore pressure, pw, at length L

if nargin < 4
    mode = 'length';
end

if strcmp(mode,'length')
    % KGD wellbore pressure
    E  = Material.solid.constitutive.constant(1);   % elastic Modulus [Pa]
    nu = Material.solid.constitutive.constant(2);   % Poisson's ratio of the solid [-]
    mu = Material.fluid.mu;                         % fluid viscosity [Pa.s]

    % In situ stress field 
    % S0  = [Domain.InsituStress.Sx,  Domain.InsituStress.Sxy;
    %        Domain.InsituStress.Sxy, Domain.InsituStress.Sy];

    % ncr = CMesh.surface_normal(1,:);   
    % S1  = S0(1,1)/2+S0(2,2)/2 + abs(S0(1,1)/2-S0(2,2)/2)*(ncr(2)^2 - ncr(1)^2);
    % S2  = S0(1,1)/2+S0(2,2)/2 - abs(S0(1,1)/2-S0(2,2)/2)*(ncr(2)^2 - ncr(1)^2);

    % Smin = min(abs(S1),abs(S2));
    % Smin = S0(1,1)/2+S0(2,2)/2 + abs(S0(1,1)/2-S0(2,2)/2)*(ncr(2)^2 - ncr(1)^2);

    q = 1e-4;       % injection flux [m^2/s]
    h = 1;          % fracture height [m]
    Q = q*h;        % injection flow rate [m^3/s]
    G = E/(1+nu)/2; % shear modulus [Pa]

    Smin = 0;
    pw = -Smin + .96*( 2*Q*mu*G^3 / (var^2*(1-nu)^3) )^(1/4);  % wellbore pressure

    C = .48*( 8*G*Q^3 / ((1-nu)*mu) )^(1/6);                   % time
    t = ( var / C ).^(3/2);
    
elseif strcmp(mode,'time')
    pmax = Control.Wellbore.Pmax;   % max pressure
    Tmax = Control.Wellbore.Tmax;   % period of loading
    if var < Tmax
        pw = pmax*sin(pi/Tmax*var); % wellbore pressure
    else
        pw = 0;
    end
    t = var;
    Q = 0;
end
    
end

