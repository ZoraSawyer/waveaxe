function conn=Tricheck(node,conn,verbose)   % Edited November 3, 2011

% FUNCTION
%
%   conn=Tricheck(node,conn,verbose)
%
% This function check wether a triangle has a negative Jacobian, and if
% so reorders it so that the the Jacobian is positive.

if ( nargin==2 )
  verbose=0;
end

if ( size(node,2)==3 )
  node=node(:,1:2);
end

count=0;
Ojac = [];                                                                 % Line 19 added to the code.
for e=1:size(conn,1)
  
  sctr=conn(e,:);
  [N,dNdxi]=LagrangeBasis('T3',[1/3 1/3]);
  detJ=det(node(sctr,:)'*dNdxi);
  if ( abs(detJ) < 0.0001 )                                                % Line 25 & 28: condition for "if" and "elseif" changed.
    Ojac = [Ojac; e];                                                      % Line 21 added to the code.
    %disp(['ZERO JACOBIAN IN ELEMENT ',num2str(e),' REMOVED'])
  elseif ( detJ < 0 )
    %disp(['NEGATIVE JACOBIAN IN ELEMENT ',num2str(e)])
    conn(e,:)=fliplr(sctr);
    count=count+1;
  end
end
conn(Ojac,:)=[];                                                           % Line 34 added to the code in order to 
                                                                           % ommite zero jacobian triangles.  
%if ( verbose )
  %disp(['TRICHECK FOUND ',num2str(count),' NEGATIVE JACOBIANS, ALL FIXED'])
end