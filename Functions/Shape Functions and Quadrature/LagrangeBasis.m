function [Nv,dNdxi]=LagrangeBasis(type,coord,dim) 
%LAGRANGEBASIS returns the lagrange interpolant basis and its gradients w.r.t the
% parent coordinate system.
%
%         [N(xi),dNdxi(xi)]=LagrangeBasis(type-order,coord,dim)
%
%   type is the toplogical class of finite element it is in the general
%   form 'topology-#of nodes' ie a three node triangel is T3 a four
%   node quadralateral is Q4 a 4 node tetrahedra is H4 a 27 node brick
%   is B27 etc
%  
%   coord is the parent coordinates at which the basis and its
%   gradients are to be evaluated at.
%
%   presently defined are L2, L3, T3, T4(cubic bubble), T6, Q4, Q9,
%   H4, H10, B8 and B27  
%
%   If dim is set to 2 then the vector representation of the N
%   matrix is returned.
%
% written by Jack Chessa
%            j-chessa@northwestern.edu
% Department of Mechanical Engineering 
% Northwestern University    
    
    if ( nargin == 2 )
        dim=1;
    end
  
    switch type
    case 'L2'  
    %%%%%%%%%%%%%%%%%%%%% L2 TWO NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2
    %
      if size(coord,2) < 1
          disp('ERROR in LagrangeBasis: coordinate needed for the L2 element')
      else
          xi=coord(1);
          N=([1-xi,1+xi]/2)';
          dNdxi=[-1;1]/2;
      end
  
    case 'L3' 
    %%%%%%%%%%%%%%%%%%% L3 THREE NODE LINE ELEMENT %%%%%%%%%%%%%%%%%%%%% 
    %  
    %    1---------2----------3
    %
    if size(coord,2) < 1
      disp('ERROR in LagrangeBasis: two coordinates needed for the L3 element')
    else
     xi=coord(1);
     N=[(1-xi)*xi/(-2);(1+xi)*xi/2;1-xi^2];
     dNdxi=[xi-.5;xi+.5;-2*xi];
    end
    
   case 'T3'
    %%%%%%%%%%%%%%%% T3 THREE NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('ERROR in LagrangeBasis: two coordinates needed for the T3 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta;xi;eta];
      dNdxi=[-1,-1;1,0;0,1];
    end
            
   case 'T4'
    %%%%%%%%%% T4 FOUR NODE TRIANGULAR CUBIC BUBBLE ELEMENT %%%%%%%%%%%% 
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         /          \
    %        /      4     \
    %       /              \
    %      /                \
    %     /                  \
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('ERROR in LagrangeBasis: two coordinates needed for the T4 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-xi-eta-3*xi*eta;xi*(1-3*eta);eta*(1-3*xi);9*xi*eta];
      dNdxi=[-1-3*eta,-1-3*xi;
	     1-3*eta, -3*xi;
	     -3*eta,   1-3*xi;
	     9*eta,   9*xi ];
    end
    
   case 'T6'
    %%%%%%%%%%%%%%%%%% T6 SIX NODE TRIANGULAR ELEMENT %%%%%%%%%%%%%%%%%%
    %   
    %               3
    %             /  \
    %            /    \
    %           /      \
    %          /        \
    %         6          5
    %        /            \
    %       /              \
    %      /                \
    %     /                  \
    %    1---------4----------2
    %
    if size(coord,2) < 2
      disp('ERROR in LagrangeBasis: two coordinates needed for the T6 element')
    else
      xi=coord(1); eta=coord(2);
      N=[1-3*(xi+eta)+4*xi*eta+2*(xi^2+eta^2);
                                  xi*(2*xi-1);
                                eta*(2*eta-1);
                              4*xi*(1-xi-eta);
                                     4*xi*eta;
                              4*eta*(1-xi-eta)];
        
      dNdxi=[4*(xi+eta)-3   4*(xi+eta)-3;
                   4*xi-1              0; 
                        0        4*eta-1;
           4*(1-eta-2*xi)          -4*xi;
                    4*eta           4*xi;
                   -4*eta  4*(1-xi-2*eta)];
    end
    
    
   case 'Q4'
    %%%%%%%%%%%%%%% Q4 FOUR NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4--------------------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1--------------------2
    %
    if size(coord,2) < 2
      disp('ERROR in LagrangeBasis: two coordinates needed for the Q4 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[ (1-xi)*(1-eta);
              (1+xi)*(1-eta);
              (1+xi)*(1+eta);
              (1-xi)*(1+eta)];
      dNdxi=1/4*[-(1-eta), -(1-xi);
		         1-eta,    -(1+xi);
		         1+eta,      1+xi;
                -(1+eta),   1-xi];
    end
    
   case 'Q9'
    %%%%%%%%%%%%%%% Q9 NINE NODE QUADRILATERIAL ELEMENT %%%%%%%%%%%%%%%%
    %
    %    4---------7----------3
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    8          9         6
    %    |                    |
    %    |                    |
    %    |                    |
    %    |                    |
    %    1----------5---------2
    %
    if size(coord,2) < 2
      disp('ERROR in LagrangeBasis: two coordinates needed for the Q9 element')
    else
      xi=coord(1); eta=coord(2);
      N=1/4*[xi*eta*(xi-1)*(eta-1);
             xi*eta*(xi+1)*(eta-1);
             xi*eta*(xi+1)*(eta+1);
             xi*eta*(xi-1)*(eta+1);
            -2*eta*(xi+1)*(xi-1)*(eta-1);
            -2*xi*(xi+1)*(eta+1)*(eta-1);
            -2*eta*(xi+1)*(xi-1)*(eta+1);
            -2*xi*(xi-1)*(eta+1)*(eta-1);
             4*(xi+1)*(xi-1)*(eta+1)*(eta-1)];
      dNdxi=1/4*[eta*(2*xi-1)*(eta-1),xi*(xi-1)*(2*eta-1);
                 eta*(2*xi+1)*(eta-1),xi*(xi+1)*(2*eta-1);
                 eta*(2*xi+1)*(eta+1),xi*(xi+1)*(2*eta+1);
                 eta*(2*xi-1)*(eta+1),xi*(xi-1)*(2*eta+1);
                -4*xi*eta*(eta-1),   -2*(xi+1)*(xi-1)*(2*eta-1);
         -2*(2*xi+1)*(eta+1)*(eta-1),-4*xi*eta*(xi+1);
                -4*xi*eta*(eta+1),   -2*(xi+1)*(xi-1)*(2*eta+1);
         -2*(2*xi-1)*(eta+1)*(eta-1),-4*xi*eta*(xi-1);
                 8*xi*(eta^2-1),      8*eta*(xi^2-1)];
    end
    
   case 'H4'
    %%%%%%%%%%%%%%%% H4 FOUR NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        /    |    \ 
    %       /     |     \
    %      1 -----|------3
    %         -   2  -
    if size(coord,2) < 3
      disp('ERROR in LagrangeBasis: three coordinates needed for the H4 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      N=[1-xi-eta-zeta;
                    xi;
                   eta;
                  zeta];
      dNdxi=[-1  -1  -1;
              1   0   0;
              0   1   0;
              0   0   1];
    end
    
   case 'H10'
    %%%%%%%%%%%%%%%% H10 TEN NODE TETRAHEDRAL ELEMENT %%%%%%%%%%%%%%%%%%
    %
    %             4
    %           / | \
    %          /  |  \
    %         /   |   \ 
    %        7    |    9 
    %       /     10     \
    %      /      |      \
    %     /       |       \
    %    1- - - -6|- - - --3
    %      -      |      -
    %        5    |    8
    %          -  |  -
    %             2
    if size(coord,2) < 3
      disp('ERROR in LagrangeBasis: three coordinates needed for the H10 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      phi=[1-xi-eta-zeta;  xi;  eta;  zeta];
      N=[phi(1)*(2*phi(1)-1);
        phi(2)*(2*phi(2)-1);
        phi(3)*(2*phi(3)-1);
        phi(4)*(2*phi(4)-1);
        4*phi(1)*phi(2);
        4*phi(1)*phi(3);
        4*phi(1)*phi(4);
        4*phi(2)*phi(3);
        4*phi(3)*phi(4);
        4*phi(2)*phi(4)];
      dNdxi = 4*[-phi(1)+.25,   -phi(1)+.25,   -phi(1)+.25;
        phi(2)-.25,    0,             0;
        0,             phi(3)-.25,    0;
        0,             0,             phi(4)-.25;
        phi(1)-phi(2), -phi(2),       -phi(2);
        -phi(3),       phi(1)-phi(3), -phi(3);
        -phi(4),       -phi(4),       phi(1)-phi(4);
        phi(3),        phi(2),        0;
        0,             phi(4),        phi(3);
        phi(4),        0,             phi(2) ];
    end
    
   case 'B8'
    %%%%%%%%%%%%%%%%%%% B8 EIGHT NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%%%%
    % 
    %                  8 
    %               /    \    
    %            /          \
    %         /                \
    %      5                     \
    %      |\                     7
    %      |   \                / |
    %      |     \     4    /     |
    %      |        \    /        |
    %      |           6          |
    %      1           |          |
    %       \          |          3
    %          \       |        /
    %            \     |     /
    %               \  |  /
    %                  2
    %                
    if size(coord,2) < 3
      disp('ERROR in LagrangeBasis: three coordinates needed for the B8 element')
    else
      xi=coord(1); eta=coord(2); zeta=coord(3);
      I1=1/2-coord/2;
      I2=1/2+coord/2;
      N=[   I1(1)*I1(2)*I1(3);
            I2(1)*I1(2)*I1(3);
            I2(1)*I2(2)*I1(3);
            I1(1)*I2(2)*I1(3);
            I1(1)*I1(2)*I2(3);
            I2(1)*I1(2)*I2(3);
            I2(1)*I2(2)*I2(3);
            I1(1)*I2(2)*I2(3)   ];
      dNdxi=[   -1+eta+zeta-eta*zeta   -1+xi+zeta-xi*zeta  -1+xi+eta-xi*eta;
                 1-eta-zeta+eta*zeta   -1-xi+zeta+xi*zeta  -1-xi+eta+xi*eta;
                 1+eta-zeta-eta*zeta    1+xi-zeta-xi*zeta  -1-xi-eta-xi*eta;
                -1-eta+zeta+eta*zeta    1-xi-zeta+xi*zeta  -1+xi-eta+xi*eta;      
                -1+eta-zeta+eta*zeta   -1+xi-zeta+xi*zeta   1-xi-eta+xi*eta;
                 1-eta+zeta-eta*zeta   -1-xi-zeta-xi*zeta   1+xi-eta-xi*eta;
                 1+eta+zeta+eta*zeta    1+xi+zeta+xi*zeta   1+xi+eta+xi*eta;
                -1-eta-zeta-eta*zeta    1-xi+zeta-xi*zeta   1-xi+eta-xi*eta  ]/8;
    end
    
   case 'B27'
    %%%%%%%%%%%%%% B27 TWENTY SEVEN NODE BRICK ELEMENT %%%%%%%%%%%%%%%%%
    %
    %                  19 
    %               /    \    
    %            20         \
    %         /                22
    %      21                    \
    %      |\         23          25
    %      |   \                / |
    %     12     24         26    |
    %      |        \    /        16
    %      |          27          |
    %      3    15     |     17   |
    %       \          |          7
    %          \       18       /
    %           6      |     8
    %               \  |  /
    %                  9
    %                
    if size(coord,2) < 3
        disp('ERROR in LagrangeBasis: three coordinates needed for the B27 element')
    else
        N=zeros(27,1);
        dNdxi=zeros(27,3);
        xi=coord(1); eta=coord(2); zeta=coord(3);
        c = 1;
        for i=1:3
            [Ni,dNdxI] = L3at(zeta,i);
            for j=1:3
                [Nj,dNdxj] = L3at(eta,j);
                for k=1:3
                    [Nk,dNdxk] = L3at(xi,k);
                    N(c) = Ni*Nj*Nk;
                    dNdxi(c,1) = Ni*Nj*dNdxk;
                    dNdxi(c,2) = Ni*dNdxj*Nk;
                    dNdxi(c,3) = dNdxI*Nj*Nk;
                    c=c+1;
                end
            end
        end       
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
    disp(['ERROR in LagrangeBasis: Element ',type,' not yet supported'])
    N=[]; dNdxi=[];
    end
 
  I = eye(dim);
  nn = size(N, 1);
  Nv = zeros(dim*nn, dim);
  for i = 1:nn
    Nv((i-1)*dim+1:i*dim, :) = I*N(i);
  end

end

function [Ni,dNdxi] = L3at(xi,index)

  N = [(1-xi)*xi/(-2);1-xi^2;(1+xi)*xi/2];
  dNdx = [xi-.5;-2*xi;xi+.5];
  Ni = N(index);
  dNdxi = dNdx(index);

end