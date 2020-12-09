function sctr = GetScatter( enodes )

global SMesh

nsd    = size(SMesh.nodes,2);     % number of space dimensions
nen    = length(enodes);          % number of enodes
entype = SMesh.EnrType(enodes);   % enrichment types of element nodes
enrH   = length(find(entype==1)); % number of nodes enriched by Heaviside function
ndof   = (nen+enrH)*nsd;          % total number of element DOFs

sctr = zeros(1,ndof);             % scatter vector
count = 0;

for n = 1:nen
    
    rangeSTD = (n-1)*nsd+1 : n*nsd;     % for 2D, 2n-1:2n  ;  for 3D, 3n-2:3n  ;  for 1D, n
    
    if (entype(n) == 1)                 % Enriched node
        count = count + 1;
        rangeENR = nen*nsd + (((count-1)*nsd+1) : count*nsd);
    else
        rangeENR = [];
    end
       
    range = [rangeSTD rangeENR];
    sctr(range) = NodalDOFs(enodes(n));  % scatter vector for the n-th enode 
end                                      
end

