function [prop_dir] = PropCriterion(inputfield,gplot)

% Compares the stress state of the domain against the fracture propagation
% criterion and returns the propagation direction for each crack tip.
%
%   Input
%           inputfield : This input argument can be either the nodal displacement
%                        field, d (number of space dimensions X 1 vector), or the nodal stress 
%                        field (number of nodes X 3 matrix). The code discreminates between these 
%                        two types of input by the field size, and uses the former to compute the
%                        actual stress field and the latter to estimate a smoothed (mapped) stress
%                        field.
%
%           gplot :      if gplot = 1, stress sampling points will be plotted 
%                        if gplot = 0, stress sampling points will not be plotted
%
%   Output
%           prop_direction : returns a ncrack by 2 matrix where ncrack is
%                            the number of crack tips
%           
%                        |x1 y1| ---> propagation direction of tip #1
%                        |x2 y2| ---> propagation direction of tip #2
%       prop_direction = | ... |
%                        |xn yn| ---> propagation direction of tip #n

% Written by Matin Parchei Esfahani, University of Waterloo, Sep. 2015


global SMesh CMesh Material

if nargin < 2
    gplot = 0;
end

ncrack = size(CMesh,2);     % number of cracks
prop_dir = zeros(ncrack,2); % propagation direction

for nc = 1:ncrack
    tip  = CMesh(nc).nodes(CMesh(nc).tip_nodes,:);  % location of the crack tip
    etip = SMesh.cmesh_e(CMesh(nc).tip_smesh_e);    % crack element containing crack tip
    n    = CMesh(nc).surface_normal(etip,:);        % unit normal vector to the surface of the last crack segment [-sin(angle) cos(angle)]
    t    = [n(2); -n(1)];                           % unit tangential vector to the surface of the last crack segment [cos(angle) sin(angle)]

    avgsize = sqrt(SMesh.eSize(CMesh(nc).tip_smesh_e));   % average element size
    %(SMesh.eSize(1)+SMesh.eSize(2))/2;    

    % Computing stress ahead of the fracture tip
    tipSpoint = tip + .001*t';    % compute stress at this point

    [stress, sp_elem] = ComputeStressatPoint(tipSpoint,inputfield);   % stress tensor at stress point

    % Calculate principal stresses at the tip
    % stress
    a = stress(1)/2 + stress(2)/2;
    b = sqrt( (stress(1)/2 - stress(2)/2)^2 + stress(3)^2 );
    % a+b
    % a-b
    Smax = max(a+b,a-b);    % maximum principal stress at tip

    if ~SMesh.einc(sp_elem) %element belongs to the domain
        T = Material.solid.TensileStrength;                                    % Tensile strength

    else                    % element belongs to an inclusion
        T  = Material.solid.inclusion(SMesh.einc(sp_elem)).TensileStrength;    % Tensile strength
    end

    % plotting stress values
    if gplot
        % plot stress values (cont'd within the loop on stress points)
        figure
        hold on
    end

    C_reduction = 0.95;
    % C_reduction*T

    if Smax > C_reduction*T     % fracture must propagate

        r     = 2.5*avgsize;                            % stress sampling radius
        theta = -71/180*pi : pi/180 : 71/180*pi;        % theta coordinate of stress points in fracture tip polar coordinates

        % Location of stress points (global)
        S = n(2).*sin(theta) - n(1).*cos(theta);
        C = n(2).*cos(theta) + n(1).*sin(theta);

        Spoints(1,:) = tip(1) + r.*C;
        Spoints(2,:) = tip(2) + r.*S;

        Stt = zeros(length(theta),1);   % hoop stress at stress points

        for npt = 1:length(theta)       % loop on stress points

            S = n(2)*sin(theta(npt)) - n(1)*cos(theta(npt));    % sin of the rotation angle
            C = n(2)*cos(theta(npt)) + n(1)*sin(theta(npt));    % cos of the rotation angle

            [stress,~] = ComputeStressatPoint(Spoints(:,npt)',inputfield);   % stress tensor at stress point and 
                                                                                % the corresponding element
            Stt(npt) = stress(1)*S^2 + stress(2)*C^2 - 2*stress(3)*C*S;         % hoop stress at stress point
            %Srt = (stress(2)-stress(1))*S*C + stress(3)*(C^2-S^2);             % shear stress at stress points

            % plot stress values (cont'd)
            if gplot
                plot(theta(npt)*180/pi,Stt(npt),'ko')
            end

        end
        if gplot
            hold off    % plot stress values
        end

        [~,sp_prop] = max(Stt); % direction of propagation

        prop_dir(nc,:) = Spoints(:,sp_prop)' - tip         % propagation direction vector (ATTENTION: Not Normalized)
    %    prop_dir = [sqrt(2)/2,sqrt(2)/2];
    %    plot stress points

        tip_elem = CMesh(nc).tip_smesh_e;       % element containing fracture tip
        tip_edge = CMesh(nc).smesh_tipedge;     % edge containing fracture tip

        switch SMesh.type
            case 'Q4'
                edge = [1,2; 2,3; 3,4; 4,1];
            case 'Q9'
                edge = [1,2; 2,3; 3,4; 4,1];
        end

        edge_nodes = SMesh.conn(tip_elem,edge(tip_edge,:));                         % nodes of the edge containing fracture tip

        edge_tangent = SMesh.nodes(edge_nodes(2),:) - SMesh.nodes(edge_nodes(1),:); % tangential vector to the edge containing fracture tip
        edge_tangent = edge_tangent/norm(edge_tangent);                             % unit tangent vector

        edge_normal = [-edge_tangent(2), edge_tangent(1)];                          % unit normal vector to the edge containing fracture tip


        
        A = dot(edge_normal,t);
        B = dot(edge_normal,prop_dir(nc,:)./norm(prop_dir(nc,:)));

        if (A*B <= 0) || (abs(B) < sin(pi/12))
            disp('Modifying fracture path to avoid instability')
%             dot(edge_tangent,t)
            sIgN1 = sign(dot(edge_tangent,t));
            sIgN2 = sign(dot(edge_normal,t)*dot(edge_tangent,t));
            prop_dir(nc,:) = sIgN1*(edge_tangent + sIgN2*tan(pi/12)*edge_normal) % Propagating with 15 degrees outwards

        end

        if gplot
            figure
            plot(tip(1),tip(2),'ro')
            axis equal
            hold on
            plot(Spoints(1,:),Spoints(2,:),'k*')
            plot(SMesh.nodes(:,1),SMesh.nodes(:,2),'bs')
            hold off
        end

    else    % Fracture does not propagate

        prop_dir(nc,:) = [0, 0];   
    end
end

end

