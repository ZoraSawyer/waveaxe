function [Mesh, Domain, Material, Control] = ConfigFile()

% INPUT/OUTPUT PATHS
    % Output files will be saved to this path (needs a "\" at the end)
    Control.OutPath = 'C:\Users\e2rivas\Documents\Matlab Results\Waveaxe\'; 
    
    % Input mesh file e.g., 1WB9in_Q4S_v2.msh must be in this folder
    Control.InPath = 'C:\Users\e2rivas\Documents\waveaxe\Config Files\';
    
    Mesh.Input = 'Gmsh';             % 'Built-in' : use the biult-in mesh generator
                                    % 'Gmsh'     : import Gmsh mesh file
                                    
% GEOMETRY AND MESH 
    Mesh.nsd = 2;                            % number of space dimensions

    if strcmp(Mesh.Input, 'Built-in')
        
        Mesh.elemtype = 'Q4';                % order of approximation. If Q4 : 4-node biliniear element
                                        %                         If Q9 : 9-node quadratic element
                                        
        Mesh.Form = 'STRUCTURED';        % STRUCTURED mesh
        
        % DOMAIN SIZE
        Domain.Lx  = 10;                       % size of domain in x-direction [m]
        Domain.Ly  = 10;                       % size of domain in y-direction [m]
        
        % MESH PARAMETERS
        Mesh.type = 'XUNIFORM';          % if "UNIFORM", uniform mesh in x- and y-directions.
                                        % if "NONUNIFORM", nonuniform mesh with ratio "rx" in
                                        % x-direction, and "ry" in y-direction.
                                        % if "XUNIFORM", uniform mesh in x-direction and nonuniform 
                                        % mesh with ratio "ry" in y-direction.
                                        % if "YUNIFORM", uniform mesh in y-direction and nonuniform 
                                        % mesh with ratio "rx" in x-direction.

        Mesh.nex = 200;                      % number of elements in x-direction (UNIFORM and XUNIFORM mesh)
        Mesh.ney = 199;                      % number of elements in y-direction (UNIFORM mesh)

        Mesh.s0x = .05;                      % finest mesh size in x-direction (NONUNIFORM mesh)
        Mesh.s0y = .05;                      % finest mesh size in y-direction (NONUNIFORM and XUNIFORM mesh)

        Mesh.rx  = 1.02;                     % mesh ratio in x-direction (NONUNIFORM mesh)
        Mesh.ry  = 1.05;                     % mesh ratio in y-direction (NONUNIFORM and XUNIFORM mesh)
        
    else
        % MESH INPUT FILE
        Mesh.FileName = '1WB9in_Q4S_v2'; % mesh file name
        Mesh.Form     = 'UNSTRUCTURED';  % STRUCTURED/UNSTRUCTURED Mesh     
    end
                                
% DOMAIN
    % MODEL CONTROLS
    Domain.Fracture_ON  = 1;            % Fracture exists  (1:Yes; 0:No)
    Domain.Inclusion_ON = 0;            % Inclusion exists (1:Yes; 0:No)
    Domain.Wellbore_ON  = 1;            % Wellbore exists  (1:Yes; 0:No)
    Domain.Leakoff_ON   = 1;            % Leak-off model   (1:Yes; 0:No)

    % IN SITU STRESS FIELD
    Domain.InsituStress.Sx  = -55.16E6; % In situ stress in x-direction [pa] (compression is negative)
    Domain.InsituStress.Sy  = -42E6; % In situ stress in y-direction [pa] (compression is negative)
    Domain.InsituStress.Sxy = 0;        % In situ shear stress [pa] (zero on principal planes)

    % OTHER DOMAIN CHARACTERISTICS
    Domain.GravityAcceleration = 9.81;      % Acceleration of Gravity [N/Kg]
    Domain.Depth               = 2459;% equivalent depth for institu pore pressure.%0.5*(10382+10736)*0.3048;%2395.88;   % Average depth [m]

    % Location of Inclusions
    if Domain.Inclusion_ON
        % Elliptical Inclusion #1
        Domain.inclusion(1).center = [.65*Lx, .50*Ly];                  % Center of the ellipse [m]
        Domain.inclusion(1).radius = [.15*Lx, .05*Ly];                  % [Rx,Ry], radii of the ellipse in x- and y-direction [m]
    end

    % WELLBORE
    if Domain.Wellbore_ON
        % Wellbore #1
        Domain.WB(1).center = [10, 10];                           % center of the wellbore
        Domain.WB(1).radiusmesh = .1222375; % radius of well in original mesh file
        Domain.WB(1).radius = 0.096837448;  % radius of the wellbore of actual.end
    end

    % FRACTURE
    if Domain.Fracture_ON
        Domain.ncrack = 2;             % number of fractures
        
        Domain.x0 = zeros(1,Domain.ncrack);   % x-coordinate of the first tip
        Domain.y0 = zeros(1,Domain.ncrack);   % y-coordinate of the first tip
        Domain.xs = zeros(1,Domain.ncrack);   % x-coordinate of the second tip
        Domain.ys = zeros(1,Domain.ncrack);   % x-coordinate of the second tip
        
        Domain.crack_surface_normal = zeros(2,Domain.ncrack);                     % normal to the negative side
        Domain.crack_front_normal   = zeros(2,Domain.ncrack);                     % normal to the first front (tip)
        
        % crack #1
        L     = .039;                                               % initial length of the fracture (must be more that one element)
        theta = 0/180*pi;                                           % inclination angle
        
        Domain.xs(1) = Domain.WB(1).center(1) + ...
                Domain.WB(1).radius*cos(theta);                     % crack start point (x-coordinate)
        Domain.ys(1) = Domain.WB(1).center(2) +... 
                Domain.WB(1).radius*sin(theta);                     % crack start point (y-coordinate)
        
        Domain.x0(1) = Domain.xs(1) + L*cos(theta);                               % location of crack tip (x-coordinate)
        Domain.y0(1) = Domain.ys(1) + L*sin(theta);                               % location of crack tip (y-coordinate)
               
        Domain.crack_surface_normal(:,1) = [-sin(theta); cos(theta)];      % normal to the negative side
        Domain.crack_front_normal(:,1)   = [cos(theta); sin(theta)];       % normal to the first front (tip)
        
        % crack #2
        L     = .039;                                               % initial length of the fracture (must be more that one element)
        theta = 180/180*pi;                                         % inclination angle
        
        Domain.xs(2) = Domain.WB(1).center(1) + ...
                Domain.WB(1).radius*cos(theta);                     % crack start point (x-coordinate)
        Domain.ys(2) = Domain.WB(1).center(2) + ...
                Domain.WB(1).radius*sin(theta);                     % crack start point (y-coordinate)
        
        Domain.x0(2) = Domain.xs(2) + L*cos(theta);                               % location of crack tip (x-coordinate)
        Domain.y0(2) = Domain.ys(2) + L*sin(theta);                               % location of crack tip (y-coordinate)
               
        Domain.crack_surface_normal(:,2) = [-sin(theta); cos(theta)];      % normal to the negative side
        Domain.crack_front_normal(:,2)   = [cos(theta); sin(theta)];       % normal to the first front (tip)
            
    else
        Domain.ncrack = 0;
        Domain.crack_surface_normal = [];
        Domain.crack_front_normal   = [];
    end

% MATERIAL PROPERTIES AND CONSTITUTIVE LAW
    % Solid
    Material.solid.rho                      = 2300;                     % density of the solid+water [kg/m^3]
    Material.solid.TensileStrength          = 300*6894.76;% 1.2E6;                      % tensile strength of the solid [pa] 
    K = 3.4e6*6894.76; % Bulk Modulus in Pa. 1 psi = 6894.76 Pa
    G = 1.4e6*6894.76; % Shear Modulus in Pa. 1 psi = 6894.76 Pa
    Material.solid.constitutive.constant(1) = 9*K*G/(3*K+G); %20E9;                     % E, Young's modulus [Pa], See DOE doc on properties of shale
    Material.solid.constitutive.constant(2) = (3*K-2*G)/2/(3*K+G); %0.20;                     % nu, Poison's ratio, See DOE doc on properties of shale
    Material.solid.constitutive.state       = 'PLANESTRAIN';            % solid stress-strain state

    % Solid - Inclusions
    if Domain.Inclusion_ON
        % Inclusion #1
        Material.solid.inclusion(1).rho                      = 2300;    % density of the solid inclusion [kg/m^3]
        Material.solid.inclusion(1).TensileStrength          = 2E6;     % tensile strength of the solid inclusion [pa] 
        Material.solid.constitutive.inclusion(1).constant(1) = 20E9;    % E, Young's modulus of the solid inclusion [Pa]
        Material.solid.constitutive.inclusion(1).constant(2) = 0.20;    % nu, poison's ratio of the solid inclusion [-]
    end

    % Solid - Cohesive law
    Material.solid.cohesive.rule  = 'BILINEAR';                         % Cohesive law
    Material.solid.FractureEnergy = 0.5*(51 + 129); %91;                                 % Gc, Fracture aperture [J/m^2]
    Material.solid.wratio         = .0005;                                 % ratio of ww/wc (BILINEAR constitutive model)

    % Fluid
    Material.fluid.rho = 1E3;                                           % density of the fluid [kg/m^3]  
    Material.fluid.mu  = 1E-3;                                          % mu, dynamic viscosity of water [Pa-s or 10^3 cP].
    Material.fluid.constitutive.wmin = 1E-3;                            % minimum aperture for fluid constitutive model near the fracture tip

    % Leak-off
    b = 7.70; % sqrt( porosity x permeability )
    phi = 0.25;
    Material.solid.porosity        = phi; % .1322;                             % porosity of the solid matrix [%]
    Material.solid.permeability    = (b^2/phi)/1000*0.9869E-12  ;% 2.56600058E-16;                    % permeability of the solid matrix [m^2]
    Material.solid.compressibility = 7.2E-11;                           % compressibility of the porous matrix [Pa^-1]

    Material.fluid.compressibility = 4.2E-10;                           % compressibility of the fluid [Pa^-1]

% COMPUTATION CONTROLS
    Control.Dynamic_ON     = 1;                                         % Dynamic analysis when 1, quasi-static analysis when 0
    Control.Wellbore.Q_inj = 1E-4;                                      % injection rate [m^3/s]
    Control.Wellbore.rate_control = 0;
    
    % pressure pulse parameters
    Control.Wellbore.npulse = 4;                                        % number of pressure pulses
    Control.Wellbore.Pmax   = 30E6;                                     % maximum pressure [Pa]
    Control.Wellbore.Tmax   = 1E-3;                                     % period of pressure pulse [s]

    Control.Time.dtmax = 4E-2;                                          % maximum time increment [s]
    Control.Time.dtmin = 1E-5;                                          % minimum time increment [s]
    Control.Time.tend  = 3E-5; % 2E-3;                                          % length of analysis [s]

    Control.coupling.maxiter   = 200;                                   % maximum number of iterations for the coupled system 
    Control.coupling.tolerance = 1E-5;                                  % admissible error of the coupled system iteration
    Control.solid.maxiter      = 300;                                   % maximum number of iterations for Newton-Raphson method
    Control.solid.tolerance    = 1E-6;                                  % admissible error of the NR process
    Control.fluid.maxiter      = 0;                                     % maximum number of iterations for Newton-Raphson method
    Control.fluid.tolerance    = 0;                                     % admissible error of the NR process

    Control.MagCoef = 1;                                                % Magnification coefficient                                    

    Control.Postprocessing.OutputFreq = 10;                              % output frequency; frequency of saving output files
    Control.Postprocessing.PlotFreq   = 10;                              % frequency of fracture plots 
    
end