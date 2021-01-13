function [M, Kuu, Kcoh, Kup, Kpu, Kpp, Fcoh, Fp, Kpp_L, FL, S11, S12] = ...
    ComputeCoupledMatrices(d, p, t, update, SMesh, CMesh, Material, Control, Domain)
% COMPUTECOUPLEDMATRICES Computes tangential matrices required for a fully coupled HF model
%   Input - d: nodal values of solid displacements
%           p: nodal values of fluid pressure
%           update: If update = 1, domain integrals will be updated. 
%                   Otherwise (update = 0), only boundary integrals will be updated.

%   Output - Kuu   = dF_int/du
%            Kcoh  = dFcoh/du
%            Kup   = dFp/dp   
%            Kpu   = dFq/du
%            Kpp   = dFq/dp
%            Kpp_L = dFL/dp
%            Fcoh: Cohesive force vector
%            Fp  : Force on the crack boundary due to fluid pressure
%            FL  : Fluid leak-off flux vector 

% Written by Matin Parchei Esfahani, May 2017, University of Waterloo
% Revised, November 2018

%in situ stress field
Sx  = Domain.InsituStress.Sx;   % in situ stress in x-direction
Sy  = Domain.InsituStress.Sy;   % in situ stress in y-direction
Sxy = Domain.InsituStress.Sxy;  % in situ shear stress (zero on principal planes)

S0 = [Sx, Sxy; Sxy, Sy];        % in situ stress tensor

nsd    = size(SMesh.nodes,2);   % number of space dimensions
ncrack = size(CMesh,2);         % number of cracks

ne_domain = size(SMesh.conn,1); % number of domain elements
ne_crack  = 0;
ne_cr     = zeros(ncrack,1);
for n = 1:ncrack
    ne_cr(n) = size(CMesh(n).conn,1);   % number of crack elements in crack mesh
    ne_crack = ne_crack + ne_cr(n);     % total number of crack elements 
end

nne_s = size(SMesh.conn,2);     % number of nodes per element (domain) 

etype = SMesh.type;             % domain element type

ndof_domain = length(d);        % number of DoFs (SMesh)
ndof_crack  = length(p);        % nember of DoFs (CMesh)

if strcmp(etype,'Q4') || strcmp(etype,'Q9')
    nEnrne = 4;                 % number of enriched nodes in an enriched element
end

% Solid model parameters
rho = Material.solid.rho;                                                   % mass density of solid

% Fluid model parameters
mu    = Material.fluid.mu;                                                  % fluid viscosity
Gc    = Material.solid.FractureEnergy;                                      % fracture energy
fu    = Material.solid.TensileStrength;                                     % ultimate cohesive traction (tensile strength)
alpha = 1;                                                                  % fluid critical aperture correction coef.
wmin  = alpha*(2*Gc/fu);                                                    % wmin = alpha * wc of the cohesive model                                                                             
wmin  = max(Material.fluid.constitutive.wmin, wmin);                        % minumum fracture aperture for fluid model

k     = Material.solid.permeability;                                        % permeability of solid matrix
phi   = Material.solid.porosity;                                            % porosity of solid matrix 
ct    = Material.solid.compressibility + Material.fluid.compressibility;    % total compressibility (ct = cf + cs)                

p_hyd = Domain.GravityAcceleration * Domain.Depth * Material.fluid.rho;     % Hydro-static pressure

% DOMAIN INTEGRALS
if update     % Domain integrals must be updated
    
    % initializing counters
    count = 1;

    % defining size of vectors
    vec_size1 = (ne_domain - ne_crack) * (nne_s * nsd)^2 + ...
                 ne_crack * ((nne_s + nEnrne) * nsd)^2;    % vector size (solid dofs)

    % pre-allocating memory         
    row = zeros(vec_size1,1);
    col = zeros(vec_size1,1);
    Kuu = zeros(vec_size1,1);       % vectorized stiffness matrix of the solid
    M   = zeros(vec_size1,1);       % vectorized mass matrix of the solid

    for e = 1:ne_domain

        switch etype                % determining order of quadrature rule
            case 'Q4'
                nq = 2;
                fLSrange = 5:8;     % range of normal LSs for each element
            case 'Q9'
                nq = 3;
                fLSrange = 5:8;     % range of normal LSs for each element
        end

        if ismember(e,SMesh.EnrElements)    % Crack element (Sub-triangles must be defined)
            nq = 5;                                                     % order of quadrature rule for each sub-triangle              
            crnum = SMesh.Crnum(e);                                     % coresponding crack mesh
            [W,Q] = DiscontQ4quad(nq,SMesh.eLS(e,fLSrange,crnum));      % quadrature points and weights
            nq = length(W);                                             % total number of quadrature points
        else                                % uncracked elements (may contain enriched nodes)
            crnum = 1;                                                  % for uncracked elements LS of the first crack is passed to Bmatrix
            [W,Q] = Quadrature(nq, 'GAUSS', nsd);                       % quadrature points and weights
            nq    = length(W);                                          % total number of quadrature points
        end

        S_enodes = SMesh.conn(e,:);         % element connectivity
        sctr     = GetScatter(S_enodes, SMesh);    % element DOFs
        S_nodes  = SMesh.nodes(S_enodes,:); % element nodal coordinates

        % Compute Kuu
        e_ndof   = length(sctr);            % number of dofs per element
        Kuu_e    = zeros(e_ndof,e_ndof);    % element stiffness matrix
        M_e      = zeros(e_ndof,e_ndof);    % element mass matrix

        D = SolidConstitutive(e, Material, SMesh, Domain);           % Solid elasticity matrix
       
        for q = 1:nq    % loop on quadrature points

           xi = Q(q,:);
           Wi = W(q);
           
           Nu = Nmatrix(xi, S_nodes, S_enodes, SMesh.EnrType(S_enodes),...
               SMesh.eLS(e,fLSrange,crnum), etype, nsd);         % N matrix at quadrature point

           [Bu, Je] = Bmatrix(xi, S_nodes, S_enodes, SMesh.EnrType(S_enodes),...
               SMesh.eLS(e,fLSrange,crnum), etype, nsd);         % B matrix and Jacobian at quadrature point

           Kuu_e = Kuu_e +  Bu' * D * Bu * Wi * det(Je);    % element stiffness matrix
           M_e   = M_e   +  Nu' * rho * Nu * Wi * det(Je);  % element mass matrix
        end

        % Forming the sparse stiffness matrix
        for i = 1:e_ndof
            col(count:count+e_ndof-1) = sctr(i)*ones(e_ndof,1);    % column index
            row(count:count+e_ndof-1) = sctr;                      % row index
            count = count + e_ndof;                                % component counter
        end

        Kuu_e = reshape(Kuu_e,[numel(Kuu_e),1]);
        Kuu(count-numel(Kuu_e):count-1) = Kuu_e;
        
        M_e = reshape(M_e,[numel(M_e),1]);
        M(count-numel(M_e):count-1) = M_e;

    end

    Kuu = sparse(row,col,Kuu,ndof_domain,ndof_domain);             % solid stiffness matrix
    M   = sparse(row,col,M,ndof_domain,ndof_domain);               % solid mass matrix
    
else    % domain integrals will not be updated
    Kuu = 0;
    M   = 0;
end

% BOUNDARY (CRACK) INTEGRALS
% defining size of vectors
nne_c     = size(CMesh(1).conn,2);                        % number of nodes per element (crack)
vec_size2 = ne_crack * nne_c^2;                           % vector size (crack dofs)
vec_size3 = ne_crack * ((nne_s + nEnrne)* nsd) * nne_c;   % vector size (coupled dofs)
vec_size4 = ne_crack * ((nne_s + nEnrne)* nsd)^2;         % vector size (cohesive dofs)
vec_size5 = ne_crack * ((nne_s + nEnrne)* nsd);           % vector size (cohesive force)
vec_size6 = ne_crack * nne_c;                             % vector size (leak-off force)

% Pre-allocating memory
Kpp     = zeros(vec_size2,1);     
Kup     = zeros(vec_size3,1);      
Kpu     = zeros(vec_size3,1);      
Kcoh    = zeros(vec_size4,1);      
Fcoh    = zeros(vec_size5,1);
Fp      = zeros(vec_size5,1);
Kpp_L   = zeros(vec_size2,1);
FL      = zeros(vec_size6,1);

% row vectors of the sparse matrices
row_pp   = zeros(vec_size2,1);
row_pu   = zeros(vec_size3,1);
row_up   = zeros(vec_size3,1);
row_coh  = zeros(vec_size4,1);
row_Fcr  = zeros(vec_size5,1);
row_FL   = zeros(vec_size6,1);

% column vectors of the sparse matrices
col_pp  = zeros(vec_size2,1);
col_pu  = zeros(vec_size3,1);
col_up  = zeros(vec_size3,1);
col_coh = zeros(vec_size4,1);

% Initializing counters
count_pp   = 1;  
count_pu   = 1;
count_up   = 1;
count_coh  = 1;
count_Fcr  = 1;
count_FL   = 1;

end_sctr = 0;   

for nc = 1:ncrack
    for e = 1:ne_cr(nc)

        switch etype
            case 'Q4'
                fLSrange = 5:8;                 % range of normal LSs for each element
            case 'Q9'
                fLSrange = 5:8;                 % range of normal LSs for each element
        end
        
        t0 = CMesh(nc).t0(e);                   % initiation time of the fracture segment
        
        nq = 5;                                 % Order of quadrature
        [W,Q] = Quadrature(nq, 'GAUSS', nsd-1); % quadrature points and weights
        nq    = length(W);                      % total number of quadrature points

        C_enodes = CMesh(nc).conn(e,:);         % crack element connectivity (CMesh)
        C_nodes  = CMesh(nc).nodes(C_enodes,:); % global coordinates of the crack nodes
        C_sctr   = C_enodes + end_sctr;         % crack element DoFs (CMesh)

        S_elem   = CMesh(nc).smesh_e(e);        % coresponding mesh element
        S_enodes = SMesh.conn(S_elem,:);        % element nodes (SMesh)
        S_nodes  = SMesh.nodes(S_enodes,:);     % nodal coordinates of the element (global)
        S_sctr   = GetScatter(S_enodes, SMesh);        % cracked element DOFs (SMesh)

        s    = CMesh(nc).CrackLength(C_enodes); % crack coordinates of nodes
        n_cr = CMesh(nc).surface_normal(e,:)';  % unit normal vector of the crack surface (negative side)

        e_ndof_s = length(S_sctr);              % element DoFs (domain)
        e_ndof_c = length(C_sctr);              % element DoFs (crack)

        % pre-allocating memory to element matrices
        Kpp_e    = zeros(e_ndof_c,e_ndof_c);
        Kup_e    = zeros(e_ndof_s,e_ndof_c);
        Kpu_e    = zeros(e_ndof_c,e_ndof_s);
        Kcoh_e   = zeros(e_ndof_s,e_ndof_s);
        Fcoh_e   = zeros(e_ndof_s,1);
        Fp_e     = zeros(e_ndof_s,1);
        Kpp_L_e  = zeros(e_ndof_c,e_ndof_c);
        FL_e     = zeros(e_ndof_c,1);
        
        for q = 1:nq

            xi = Q(q,:);
            Wi = W(q); 

            Np = [.5*(1-xi) .5*(1+xi)];         % fluid N matrix at quadrature point
            Je = [-.5 .5]*s;                    % Jacobian
            Bp = Je\[-.5 .5];                   % fluid B matrix at quadrature point


            % Compute jump in shape functions, Ndis = N|pos - N|neg
            X = Np * C_nodes;                             % global coordinates of the quadrature point
            xi = ParentCoordinates(X, etype, S_nodes, SMesh.Form);    % Parent coordinates of the quadrature point

            
            N_pos = Nmatrix(xi, S_nodes, S_enodes, SMesh.EnrType(S_enodes),...
                SMesh.eLS(S_elem,fLSrange,nc), etype, nsd, 1);     % Shape functions on the positive side of the crack
            N_neg = Nmatrix(xi, S_nodes, S_enodes, SMesh.EnrType(S_enodes),...
                SMesh.eLS(S_elem,fLSrange,nc), etype, nsd, -1);    % Shape functions on the negative side of the crack

            Ndis = N_pos - N_neg;   % Jump!            

            % Compute aperture, Dw, dDw/dw, and dt_coh/dw
            w = n_cr'*Ndis*d(S_sctr);           % fracture aperture at GP
            [tcoh, dtcohdw] = GetCohesive(w, Material);    % dt_coh/dw
            w = max(w,wmin);                    % modified aperture for cubic law
            Dw = w^3/mu/12;                     % Poiseuille flow (cubic law)
            dDwdw = w^2/mu/4;                   % dDw/dw

            % Compute pressure gradient
            p_gp = Np*p(C_sctr);            % pressure at quadrature point
            
            dpds = Bp*p(C_sctr);                % pressure gradient at quadrature point
            
            % Compute Leak-off coefficient
            if Domain.Leakoff_ON && (t > t0)
                C_leakoff = 2*sqrt(k*phi*ct/(mu*pi))/sqrt(t-t0);     % Carter leak-off coefficient
            else
                C_leakoff = 0;
            end          
            
            % Compute Kpp
            Kpp_e = Kpp_e + Bp' * Dw * Bp * Wi * det(Je);
            
            % Compute Kup
            Kup_e = Kup_e + Ndis' * n_cr * Np * Wi * det(Je);

            % Compute Kpu
            Kpu_e = Kpu_e + Bp' * dpds * dDwdw * n_cr' * Ndis * Wi * det(Je);
                        
            % Compute Kcoh
            Kcoh_e = Kcoh_e + Ndis' * n_cr * dtcohdw * n_cr' * Ndis * Wi * det(Je);

            % Compute Fcoh
            Fcoh_e = Fcoh_e + Ndis' * (tcoh*eye(nsd) * n_cr) * Wi * det(Je);
            
            % Compute Fp
            Fp_e = Fp_e + Ndis' * ((p_gp*eye(nsd) + S0) * n_cr) * Wi * det(Je);
                    
            % Compute Kpp_L
            Kpp_L_e = Kpp_L_e + Np' * C_leakoff * Np * Wi * det(Je);
            
            % Compute FL
            FL_e = FL_e + Np' * C_leakoff * (p_gp - p_hyd) * Wi * det(Je);

        end

        % Forming vectorized matrices
        % Compute Kpp and Kpp_L
        for i = 1:e_ndof_c
            col_pp(count_pp:count_pp+e_ndof_c-1) = C_sctr(i)*ones(e_ndof_c,1);
            row_pp(count_pp:count_pp+e_ndof_c-1) = C_sctr;
            count_pp = count_pp + e_ndof_c;
        end

        Kpp_e = reshape(Kpp_e,[numel(Kpp_e),1]);
        Kpp(count_pp-numel(Kpp_e):count_pp-1) = Kpp_e;
        
        Kpp_L_e = reshape(Kpp_L_e,[numel(Kpp_L_e),1]);
        Kpp_L(count_pp-numel(Kpp_L_e):count_pp-1) = Kpp_L_e;
        
        % Compute Kup
        for i = 1:e_ndof_c
            col_up(count_up:count_up+e_ndof_s-1) = C_sctr(i)*ones(e_ndof_s,1);
            row_up(count_up:count_up+e_ndof_s-1) = S_sctr;
            count_up = count_up + e_ndof_s;
        end

        Kup_e = reshape(Kup_e,[numel(Kup_e),1]);
        Kup(count_up-numel(Kup_e):count_up-1) = Kup_e;

        % Compute Kpu
        for i = 1:e_ndof_s
            col_pu(count_pu:count_pu+e_ndof_c-1) = S_sctr(i)*ones(e_ndof_c,1);
            row_pu(count_pu:count_pu+e_ndof_c-1) = C_sctr;
            count_pu = count_pu + e_ndof_c;
        end

        Kpu_e = reshape(Kpu_e,[numel(Kpu_e),1]);
        Kpu(count_pu-numel(Kpu_e):count_pu-1) = Kpu_e;

        % Compute Kcoh
        for i = 1:e_ndof_s
            col_coh(count_coh:count_coh+e_ndof_s-1) = S_sctr(i)*ones(e_ndof_s,1);
            row_coh(count_coh:count_coh+e_ndof_s-1) = S_sctr;
            count_coh = count_coh + e_ndof_s;
        end

        Kcoh_e = reshape(Kcoh_e,[numel(Kcoh_e),1]);
        Kcoh(count_coh-numel(Kcoh_e):count_coh-1) = Kcoh_e;

        % Compute Fcoh and Fp
        row_Fcr(count_Fcr:count_Fcr+e_ndof_s-1) = S_sctr;
        count_Fcr = count_Fcr + e_ndof_s;

        Fcoh(count_Fcr-length(Fcoh_e):count_Fcr-1) = Fcoh_e;
        Fp(count_Fcr-length(Fcoh_e):count_Fcr-1)   = Fp_e;
        
        % Compute FL
        row_FL(count_FL:count_FL+e_ndof_c-1) = C_sctr;
        count_FL = count_FL + e_ndof_c;
        
        FL(count_FL-length(FL_e):count_FL-1) = FL_e;

    end
    end_sctr = end_sctr + size(CMesh(nc).nodes,1);     % last fluid dof of the crack
    
end

% Forming sparse matrices
Kpp     = sparse(row_pp, col_pp, Kpp, ndof_crack, ndof_crack);
Kup     = sparse(row_up, col_up, Kup, ndof_domain, ndof_crack);
Kpu     = sparse(row_pu, col_pu, Kpu, ndof_crack, ndof_domain);
Kcoh    = sparse(row_coh, col_coh, Kcoh, ndof_domain, ndof_domain);
Fcoh    = sparse(row_Fcr, ones(length(row_Fcr),1), Fcoh, ndof_domain, 1);
Fp      = sparse(row_Fcr, ones(length(row_Fcr),1), Fp, ndof_domain, 1);
Kpp_L   = sparse(row_pp, col_pp, Kpp_L, ndof_crack, ndof_crack);
FL      = sparse(row_FL, ones(length(row_FL),1), FL, ndof_crack, 1);

end

