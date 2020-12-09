function [Q, Q_avg] = ComputeFlowrate(t, dt, Q_avg_old)
% COMPUTEFLOWRATE Computes rate of change in fracture volume
%
% Input
%       t         : time
%       dt        : time increment
%       Q_avg_old : average flow rate until the previous time step
%
% Output
%       Q     : instantaneous flow rate at time t
%       Q_avg : average flow rate until time t

% Written by Matin Parchei Esfahani, University of Waterloo, July 2017


global CMesh

ncrack = size(CMesh,2);     % number of fractures
Vcr = zeros(ncrack,1);      % volume of fracture
Q_avg = zeros(ncrack,1);    % average flowrate at time t
Q = zeros(ncrack,1);        % flowrate at time t

dif = length(Q_avg) - length(Q_avg_old);
if dif > 0
    Q_avg_old = [Q_avg_old; zeros(dif,1)];
end

for nc = 1:ncrack
    ne  = size(CMesh(nc).conn,1);                                               % number of fracture elements
    glc = CMesh(nc).GLconn;                                                     % global connectivity of the crack
    for e = 1:ne                                                                % compute fracture volume (trapezoidal rule)
        ds = CMesh(nc).CrackLength(glc(e+1))-CMesh(nc).CrackLength(glc(e));     % segment length
        sw = CMesh(nc).w(glc(e)) + CMesh(nc).w(glc(e+1));
        Vcr(nc) = Vcr(nc) + sw*ds/2;                                            % fracture volume
    end
    Q_avg(nc) = Vcr(nc) / t;                                                    % average flowrate
    Q(nc)     = ( Vcr(nc) - Q_avg_old(nc)*(t-dt) ) / dt;                        % flowrate
end

end

