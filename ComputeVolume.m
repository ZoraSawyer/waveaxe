function [Vcr, Vinj] = ComputeVolume(dt, Vinj_old, Qinj)

% Computes rate of change in fracture volume
%
% Input
%       dt        : time increment
%       V_inj_old : Volume injected until the previous time step
%
% Output
%       V     : Fracture storage volume at time t
%       V_inj : Volume injected until time t

% Written by Matin Parchei Esfahani, University of Waterloo, July 2017


global CMesh

ncrack = size(CMesh,2);     % number of fractures
Vcr    = zeros(ncrack,1);   % volume of fracture
Vinj   = zeros(ncrack,1);   % average flowrate at time t

dif = length(Vinj) - length(Vinj_old);
if dif > 0
    Vinj_old = [Vinj_old; zeros(dif,1)];
end

for nc = 1:ncrack
    ne  = size(CMesh(nc).conn,1);                                               % number of fracture elements
    glc = CMesh(nc).GLconn;                                                     % global connectivity of the crack
    for e = 1:ne                                                                % compute fracture volume (trapezoidal rule)
        ds = CMesh(nc).CrackLength(glc(e+1))-CMesh(nc).CrackLength(glc(e));     % segment length
        sw = CMesh(nc).w(glc(e)) + CMesh(nc).w(glc(e+1));
        Vcr(nc) = Vcr(nc) + sw*ds/2;                                            % fracture volume
    end
    Vinj(nc) = Vinj_old(nc) + Qinj(nc)*dt;                                      % volume injected
end

end

