function [ elem ] = FindElement(X, SMesh)
% Finds the element containing point X

%   Input
%        X : global coordinates of the point of interest
%   Output
%        elem : the element containing point X

% Written by Matin Parchei Esfahani, University of Waterloo, June 2015
% Last modified Oct. 2017.

ne = size(SMesh.conn,1);    % number of elements

for e = 1:ne
    enodes = SMesh.conn(e,:);       % element connectivity
    x_elem = SMesh.nodes(enodes(1:4),1); % X-coordinate of element nodes
    y_elem = SMesh.nodes(enodes(1:4),2); % Y-coordinate of element nodes
    
    [IN,ON] = inpolygon(X(1), X(2), x_elem, y_elem);
    
    if IN || ON     % point X is located either inside the e-th element or on its boundary
        elem = e;
        break
    end
end

end

