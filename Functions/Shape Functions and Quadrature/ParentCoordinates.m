function [pt] = ParentCoordinates(Gpt,elemType,nodes)

% Calculates local coordinates of a point from its global coordinates.

global SMesh

if (strcmp(elemType,'Q4') || strcmp(elemType,'Q8') || strcmp(elemType,'Q9')) ...
        && strcmp(SMesh.MeshForm, 'STRUCTURED')
    
    Xbar = sum(nodes(1:4,:))/4;                                            % Parent coordinates (xi-eta) could be
    l = norm( nodes(2,:)-nodes(1,:) );                                     % exactly calculated for quadralateral
    w = norm( nodes(3,:)-nodes(2,:) );                                     % (Qn) elements. (added to the code)
    pt = 2*(Gpt - Xbar)./[l w];
else
    coord   = zeros(1,2);
    ksi     = 0;                                                             
    eta     = 0;
    iter    = 10;
    epsilon = 0.00001;
    inc     = 1;
    while (inc < iter)  % Newton-Raphson process for local coordinates in other element types.
        [N,dNdxi]=LagrangeBasis(elemType,coord);   % compute shape functions
                                                                           
        x = N'*nodes(:,1);                                                 
        y = N'*nodes(:,2);
        
        df1dr = dNdxi(:,1)' * nodes(:,1);
        df1ds = dNdxi(:,2)' * nodes(:,1);
        df2dr = dNdxi(:,1)' * nodes(:,2);
        df2ds = dNdxi(:,2)' * nodes(:,2);
 
        f1 = x - Gpt(1);
        f2 = y - Gpt(2);

        detF = df1dr*df2ds - df1ds*df2dr ;

        invf(1,1) =  1.0/detF * df2ds;
        invf(1,2) = -1.0/detF * df1ds;
        invf(2,1) = -1.0/detF * df2dr;
        invf(2,2) =  1.0/detF * df1dr;

        ksi = ksi - invf(1,1)*f1 - invf(1,2)*f2;
        eta = eta - invf(2,1)*f1 - invf(2,2)*f2;

        coord(1) = ksi;
        coord(2) = eta;

        if( (abs(ksi - coord(1)) < epsilon) && ...
                (abs(eta - coord(2)) < epsilon) )
            inc  = iter + 1;
            pt = coord;
        else
            inc = inc + 1;
        end
    end
end
end

