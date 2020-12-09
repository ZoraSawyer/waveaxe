% Reads crack length, location of the pysical tip and time from
% "TipLocation.dat" and fracture aperture and fluid pressure from
% "Fracture.vtk" and plot fracture length, normalized cohesive zone length,
% fracture profile and pressure profile

% Written by Matin Parchei Esfahani, University of Waterloo, Sep. 2017

global IOPath CMesh

nplot = 5;                 % number of curves for each figure
dn = floor((nt-1)/nplot);   % plot frequency

for nc = 1:ncrack
    
    filename = [IOPath 'TipLocation' num2str(nc) '.dat'];
    fileID = fopen(filename,'r');
    count = 1;
    values = zeros(nt-1,5);
    
    while ~feof(fileID)
        values(count,:) = fscanf(fileID, '%f %f %f %f %f\n', [1,5]);
        count = count + 1;
    end
    fclose(fileID);

    t      = values(:,1);   % time
    Length = values(:,4);   % crack length
    ptip   = values(:,5);   % location of the physical tip
    
    % plot tip location
    figure;
    hold on
    plot(t,Length,'b-','LineWidth',1.5);
    plot(t,ptip,'b--','LineWidth',1.5);
    hold off
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('fracture length (m)','FontSize', 20, 'FontWeight', 'bold')
    xlabel('time (s)','FontSize', 20, 'FontWeight', 'bold')
    legend('fracture tip', 'physical tip')
    
    % plot cohesive zone length
    figure;
    hold on
    plot(t,(Length-ptip)./Length,'b-','LineWidth',1.5);
    hold off
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('l_{c} / L','FontSize', 20, 'FontWeight', 'bold')
    xlabel('time (s)','FontSize', 20, 'FontWeight', 'bold')
            
    % plot fracture profile
    h1 = figure;
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('fracture aperture (mm)','FontSize', 20, 'FontWeight', 'bold')
    xlabel('fracture length (m)','FontSize', 20, 'FontWeight', 'bold')
    
    % plot pressure profile
    h2 = figure;
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('fluid pressure (MPa)','FontSize', 20, 'FontWeight', 'bold')
    xlabel('fracture length (m)','FontSize', 20, 'FontWeight', 'bold')
       
    for i = 1:dn:nt-1
        filename = [IOPath 'Fracture' num2str(nc) '.vtk.' num2str(i)];
        fileID = fopen(filename,'r');
        
        s = textscan(fileID, '%s', 'delimiter', '\n');
        p_line = find(strcmp(s{1}, 'SCALARS pressure float 1'), 1, 'first');
        w_line = find(strcmp(s{1}, 'SCALARS aperture float 1'), 1, 'first');
                
        fclose(fileID);
        
        npoint = w_line - p_line - 2;       % total number of available records for each parameter
        
        p = str2double( s{1}(p_line + 2 : 2 : p_line + 2 + npoint - 1) );
        w = str2double( s{1}(w_line + 2 : 2 : w_line + 2 + npoint - 1) );
        
        glc = CMesh.GLconn(1:npoint/2);     % global connectivity of the nodes
        
        figure(h1)
        hold on
        plot(CMesh(nc).CrackLength(glc),w(glc)*1e3,'b-','LineWidth',1)
        hold off
        
        figure(h2)
        hold on
        plot(CMesh(nc).CrackLength(glc),p(glc)/1e6,'b-','LineWidth',1)
        hold off
    end
    
end