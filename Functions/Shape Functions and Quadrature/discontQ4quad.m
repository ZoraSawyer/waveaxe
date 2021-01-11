function [W,Q]=DiscontQ4quad(order,phi) % Edited by Matin Parchei Esfahani (November 3, 2011)                                 

corner = [1 2 3 4 1];
node  = [-1 -1; 1 -1; 1 1; -1 1];
upreg=[];                                                                  % Lines 6,7,21,22 and 12 to 16 added to make crack compatible triangles.       
dnreg=[];   
% loop on element edges
for i = 1 : 4
    n1 = corner(i);
    n2 = corner(i+1);
    if ( phi(n1) > 0 )
        upreg=[upreg; node(i,:)];   % upper crack region nodes
    else
        dnreg=[dnreg; node(i,:)];   % lower crack region nodes
    end
    if ( phi(n1)*phi(n2) < 0 )
        r    = phi(n1)/(phi(n1)-phi(n2));
        pnt  = (1-r)*node(n1,:)+r*node(n2,:);
        %node = [node;pnt];
        upreg=[upreg; pnt];
        dnreg=[dnreg; pnt];
    end
end
% get decompused triangles
%node = [node;-0.75 1 ; -0.5 1 ;-0.25 1 ; 0 1 ;0.25 1 ; 0.5 1 ;0.75 1 ; -0.75 -1 ;-0.5 -1 ;-0.25 -1 ; 0 -1 ;0.25 -1 ; 0.5 -1;0.75 -1 ; -1 -0.75 ;-1 -0.5 ;-1 -0.25 ; -1 0.25;-1 0.5 ;-1 0.75 ; 1 -0.75 ;1 -0.5 ;1 -0.25 ; 1 0.25 ; 1 0.5;1 0.75 ];
%node = unique(node,'rows');                                               % Lines 27: not required
%tri = delaunay(node(:,1),node(:,2));                                      % Lines 28 to 30 Changed to lines 31 to 36 
%tri = Tricheck_new(node,tri);
utri = delaunay(upreg(:,1),upreg(:,2)); %upper crack subtriangles
utri = Tricheck(upreg,utri);
dtri = delaunay(dnreg(:,1),dnreg(:,2)); %lower crack subtriangles
dtri = Tricheck(dnreg,dtri);
tri  = [utri; dtri + size(upreg,1)]; % total subtriangles
node = [upreg; dnreg];
% loop over subtriangles to get quadrature points and weights
pt = 1;
for e = 1:size(tri,1)
    [w,q] = quadrature(order,'TRIANGULAR',2);
    % transform quadrature points into the parent element
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord = [coord(2,:);coord(1,:);coord(3,:)];
        a = det([coord,[1;1;1]])/2;
    end
    if ( a~=0 )
        for n=1:length(w)
            N = LagrangeBasis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;
            pt = pt+1;
        end
    end
end
% figure
% hold on
% triplot(tri,node(:,1),node(:,2));
% n1 = plot(Q(:,1),Q(:,2),'r*');
% % plot(cr(:,1),cr(:,2),'red')
% hold off
end








