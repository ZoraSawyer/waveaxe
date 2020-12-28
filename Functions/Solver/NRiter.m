function [Pvar,q,NR,NRs,NRf,converged] = NRiter(Force, t, dt, ndof_s, ndof_c, fixed_dofs, Pvar0, Pvar_1, dynamic_ON)
%NRITER Newton-Raphson iterative process: Solves the nonlinear coupled system of 
% equations for simulation of HF. Solves F(x)-f = 0
%   Input - Force = [Ft + Fb - F_int    % solid force vectors
%                           q       ]   % fluid flux vector
%           t         : current time
%           dt        : time increament
%           ndof_s    : Total number of DoFs for the solid domain
%           ndof_c    : Total number of DoFs for the fluid domain
%           fixed_dofs: DoFs related to essential boundary conditions (for
%                       both fluid and solid)
%           Pvar0     : Values of primary variables at previous time step, 
%                       Pvar0 = Pvar(t_n) = [u0
%                                            p0]

%   Output - Pvar : Values of primary variables at time step t_n+1, 
%                   Pvar = Pvar(t_n+1) = [u_n+1
%                                         p_n+1]
%            NR   : Norm of residual |R|
%
%            converged : returns 1 if process is completed successfully.
%                        returns 0 if convergence fails.

% Written by Matin Parchei Esfahani, May 2017, University of Waterloo

global Control

max_iter = Control.coupling.maxiter;      % Maximum iteration
res_tol  = Control.coupling.tolerance;    % Residual norm tolerance

Pvar  = Pvar0;                      % initializing primary variables
dPvar = zeros(size(Pvar));          % change of the primary variable

iter = 1;                           % iteration counter
converged = 0;                      % convergence indicator

s_dof = 1:ndof_s;                   % solid dofs
c_dof = ndof_s+1 : ndof_s+ndof_c;   % fluid dofs

free_dofs = setdiff(1:ndof_s+ndof_c, fixed_dofs);   % free dofs

R = zeros(length(Pvar0),1);         % residual of the coupled system

NR  = zeros(max_iter,1);             % norm of residual
NRs = zeros(max_iter,1);
NRf = zeros(max_iter,1);

if dynamic_ON
    disp([num2str(toc), ': Dynamic Analysis'])
else
    disp([num2str(toc), ': Quasi-static Analysis'])
end

msg_len = 1;
disp([num2str(toc), ': Newton-Raphson Iteration'])
fprintf('ITERATION  ');

update = 1;

while (iter <= max_iter) && ~converged
    
    fprintf(repmat('\b',1,msg_len))
    fprintf('%d:',iter);
    msg_len = numel(num2str(iter))+1;
    
    % Computer tangential matrices
    [Mu, Ku, Kcoh, Kup, Kpu, Kpp, Fcoh, Fp, Kpp_L, FL, Spp, Spu] = ComputeCoupledMatrices(...
            Pvar(s_dof), Pvar(c_dof), Pvar0(s_dof), Pvar0(c_dof), t, update);
    
    if update
        Kuu = Ku;
        M   = Mu;
        update = 0; % Kuu and M will not be updated in this iteration again
    end
    
    % Jacobian of the fully coupled system
    Kt = [dynamic_ON*(1/dt^2)*M + Kuu + Kcoh,     -Kup;
                             1/dt*Kup' + Kpu,      Kpp + Kpp_L];
                           
    % Compute residuals
    R(s_dof) = dynamic_ON*1/dt^2*M*(Pvar(s_dof) - 2*Pvar0(s_dof) + Pvar_1(s_dof)) + ...
        Kuu*Pvar(s_dof) + Fcoh - Fp - Force(s_dof);                                         % Ru
    R(c_dof) = 1/dt*Kup'*(Pvar(s_dof)-Pvar0(s_dof)) + Kpp*Pvar(c_dof) + Force(c_dof) + FL;  % Rp
    R(fixed_dofs) = -Kt(fixed_dofs,free_dofs)*dPvar(free_dofs);
            
    % Check convergence
    NRs(iter) = norm(R(s_dof));
    NRf(iter) = norm(R(c_dof));
    NR(iter)  = norm(R);
    
    if NR(iter) < res_tol % L2 norm of the residual
        converged = 1;       
        fprintf(' CONVERGED');
    else
        
        % Solve the system of equations
        dPvar(free_dofs) = -Kt(free_dofs,free_dofs)\R(free_dofs);
                
        % impose essential BC
        dPvar(fixed_dofs) = 0;
        
        % Update variables
        Pvar = Pvar + dPvar;
    end
    
    iter = iter + 1;
end

q = - (1/dt*Kup'*(Pvar(s_dof)-Pvar0(s_dof)) + Kpp*Pvar(c_dof));     % fluid flux

if iter > max_iter
    fprintf(' CONVERGENCE FAILED \n');
    return
end

end

