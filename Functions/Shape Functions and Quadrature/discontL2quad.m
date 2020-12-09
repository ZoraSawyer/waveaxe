function [W,Q] = DiscontL2quad(order,segment)
% Returns quadrature points for each segment of a discontinuous element
% Written by Matin Parchei Esfahani, University of Waterloo, February 2018

nsegment = length(segment)-1;   % number of segments
Q = zeros(nsegment*order,1);
W = zeros(nsegment*order,1);
count = 1;

for i=1:nsegment
    [w,q] = Quadrature(order,'GAUSS',1);    % quadrature points of the segment
    for nq = 1:length(w)
        [Nv,~] = LagrangeBasis('L2',q(nq),1);
        Q(count) = segment([i,i+1])*Nv;
        W(count) = w(nq);
        count = count + 1;
    end
end

end

