function [SMesh, CMesh, Pvar, Pvar0, Pvar_1, prop] = ...
	Propagate(SMesh, CMesh, Material, Domain)

	fdof = [CMesh.fdof];

	% ========================= COMPUTE SOLID STRESS ==========================
	% Compute Nodal Stresses
	fprintf('\t%.2f: Computing stress in solid\n', toc);
	S = ComputeNodalStress(Pvar(1:SMesh.ndof), SMesh, Material, Domain);   % Compute stress at nodes

	% in situ stress
	Sx  = Domain.InsituStress.Sx;   % in situ stress in x-direction
	Sy  = Domain.InsituStress.Sy;   % in situ stress in y-direction
	Sxy = Domain.InsituStress.Sxy;  % in situ shear stress (zero on principal planes)

	S0 = [Sx; Sy; Sxy];     % in situ stress matrix (voigt)

	S_tot = S + repmat(S0,1,size(S,2));   % total stress (S + S0)

	% ======================== FRACTURE PROPAGATION ===========================
	% Check for fracture propagation
	fprintf('\t%.2f: Evaluating fracture propagation criterion\n', toc);
	if ~isempty(CMesh(1).conn)   % at least one fracture exists
	    % CHECKING FOR FRACTURE PROPAGATION
	    gplot = 0;  % plot mode is on when gplot = 1; used for debugging purposes
	    prop_dir = PropCriterion(S_tot, gplot, SMesh, CMesh, Material);    % Direction of fracture propagation
	else                         % no fractrue exists
	    prop_dir = [0 0];
	end

	% Update parameters
	if norm(prop_dir) % fracture propagation
	    
	    % saving old values
	    s_dof_old  = SMesh.ndof;
	    fdof_old   = fdof;
	    Pvar_old   = Pvar;
	    Pvar0_old  = Pvar0;
	    Pvar_1_old = Pvar_1;
	    
	    [SMesh, CMesh] = PropagateCracks(prop_dir, n, t, SMesh, CMesh, Control.OutPath);       % Updating levelsets and enriched DOFs
	                                           % based on the new fracture configuration

	    SMesh.enrDOFs = SMesh.nsd*length(SMesh.EnrNodes);  % updating number of enriched DOFs
	    SMesh.ndof = SMesh.stdDOFs + SMesh.enrDOFs;           % updating total number of solid DOFs
	    
		f_dof  = sum([CMesh.fdof]);                     
	    
	    % Redefining vectors based on updated DOFs
	    Pvar   = zeros(SMesh.ndof+f_dof,1);              % new Pvar = [d; p1; p2; ... pn]
	    Pvar0  = zeros(SMesh.ndof+f_dof,1);              % new Pvar0
	    Pvar_1 = zeros(SMesh.ndof+f_dof,1);              % new Pvar_1
	    
	    Pvar(1:s_dof_old)   = Pvar_old(1:s_dof_old);    % d
	    Pvar0(1:s_dof_old)  = Pvar0_old(1:s_dof_old);   % d0
	    Pvar_1(1:s_dof_old) = Pvar_1_old(1:s_dof_old);  % d_1
	    
	    count     = 0;
	    count_old = 0;
	    for nc = 1:Domain.ncrack
	        Pvar(SMesh.ndof+count+1 : SMesh.ndof+count+fdof_old(nc)) = ...
	        Pvar_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));    % p
	        Pvar0(SMesh.ndof+count+1 : SMesh.ndof+count+fdof_old(nc)) = ...
	        Pvar0_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));   % p0
	        Pvar_1(SMesh.ndof+count+1 : SMesh.ndof+count+fdof_old(nc)) = ...
	        Pvar_1_old(s_dof_old+count_old+1 : s_dof_old+count_old+fdof_old(nc));  % p_1
	    
	        % padding fluid pressure in the new crack segment with the value of
	        % tip pressure to avoid pressure discontinuity
	        Pvar(SMesh.ndof+count+fdof_old(nc)+1 : SMesh.ndof+count+fdof(nc)) = ...
	            Pvar_old(s_dof_old+count_old+fdof_old(nc));
	        Pvar0(SMesh.ndof+count+fdof_old(nc)+1 : SMesh.ndof+count+fdof(nc)) = ...
	            Pvar0_old(s_dof_old+count_old+fdof_old(nc));
	        Pvar_1(SMesh.ndof+count+fdof_old(nc)+1 : SMesh.ndof+count+fdof(nc)) = ...
	            Pvar_1_old(s_dof_old+count_old+fdof_old(nc));
	        
	        count = count + CMesh(nc).fdof;
	        count_old = count_old + fdof_old(nc);
	    end

	    prop = 1;   % crack is updated
	
	else
	    prop = 0;   % crack is not updated
	end

end