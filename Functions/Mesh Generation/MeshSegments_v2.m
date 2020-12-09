function [sx,sy,nex,ney] = MeshSegments_v2(meshtype, etype, Lx, Ly, nex, ney, s0x, s0y, rx, ry)
% Returns mesh segments for a Q4 structured mesh

% Output:
%   Sx  : mesh segments in x-direction
%   Sy  : mesh segments in y-direction
%   nex : number of elements in x-direction
%   ney : number of elements in y-direction

% Written by Matin Parchei Esfahani, University of Waterloo, 2017

switch meshtype
    case 'UNIFORM'
        % x-direction
        dLx = Lx / nex;
        sx  = 0:dLx:Lx;
        
        % y-direction
        dLy = Ly / ney;
        sy  = 0:dLy:Ly;
        
    case 'NONUNIFORM'
        % x-direction
        L   = Lx;
        m   = log10(1-L/s0x*(1-rx)) / log10(rx);
        n   = 1:m;
        sx  = [0, s0x .* (1-rx.^n) ./ (1-rx)];
        sx(end) = Lx;
        nex = length(sx)-1;
        
        % y-direction
        L   = Ly/2;
        m   = log10(1-L/s0y*(1-ry)) / log10(ry);
        n   = 1:m;
        seg = [0, s0y .* (1-ry.^n) ./ (1-ry)];
        top = Ly/2 + s0y/2 + seg;
        bot = fliplr(Ly/2 - s0y/2 - seg);
        sy  = [bot, top];
        sy(1)   = 0;
        sy(end) = Ly;
        ney = length(sy)-1;
        
    case 'XUNIFORM'
        % x-direction
        dLx = Lx / nex;
        sx  = 0:dLx:Lx;
        
        % y-direction
        L   = Ly/2;
        m   = log10(1-L/s0y*(1-ry)) / log10(ry);
        n   = 1:m;
        seg = [0, s0y .* (1-ry.^n) ./ (1-ry)];
        top = Ly/2 + s0y/2 + seg;
        bot = fliplr(Ly/2 - s0y/2 - seg);
        sy  = [bot, top];
        sy(1)   = 0;
        sy(end) = Ly;
        ney = length(sy)-1;
        
    case 'YUNIFORM'
        % x-direction
        L   = Lx/2;
        m   = log10(1-L/s0x*(1-rx)) / log10(rx);
        n   = 1:m;
        seg = [0, s0x .* (1-rx.^n) ./ (1-rx)];
        top = Lx/2 + s0x/2 + seg;
        bot = fliplr(Lx/2 - s0x/2 - seg);
        sx  = [bot, top];
        sx(1)   = 0;
        sx(end) = Lx;
        nex = length(sx)-1;
        
        % y-direction
        dLy = Ly / ney;
        sy  = 0:dLy:Ly;
end

if strcmp(etype,'Q9')
    
    % adding nodes in x-direction
    s1 = sx;
    s2 = (sx(1:end-1) + sx(2:end))/2;
    sx = zeros(1,2*nex+1);
    sx(1:2:end)   = s1;
    sx(2:2:end-1) = s2;
    
    % adding nodes in y-direction
%     s1 = sy;
%     s2 = (sy(1:end-1) + sy(2:end))/2;
    s1 = [sy(1), sy(2:end-1)+s0y/4, sy(end)];
    s2 = (s1(1:end-1) + s1(2:end))/2;
    sy = zeros(1,2*ney+1);
    sy(1:2:end)   = s1;
    sy(2:2:end-1) = s2;
        
end
% V = ones(length(sx), length(sy));
% surface(sx, sy, V')
% axis equal
end

