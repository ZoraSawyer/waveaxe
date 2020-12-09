function [ dof ] = NodalDOFs( node )
global SMesh

nsd      = size(SMesh.nodes,2);   % number of space DOFs
nn       = size(SMesh.nodes,1);   % number of nodes
entype   = SMesh.EnrType(node);   % type of enrichment

dofSTD = (node-1)*nsd+1:node*nsd;  % Standard FEM DOFs

if entype == 1  % Enrichment DOFs
    nenrnode = find(SMesh.EnrNodes == node);
    dofENR = nn*nsd + (((nenrnode-1)*nsd+1):nenrnode*nsd);
else
    dofENR = [];
end

dof = [dofSTD dofENR];  % Total DOFs of the node

end