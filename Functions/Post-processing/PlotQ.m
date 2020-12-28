% Reads data from the output file "Flowrate" and plots instantaneous and
% average flow rates

% Written by Matin Parchei Esfahani, University of Waterloo, July 2017

for nc = 1:ncrack
    filename = [OutPath 'Flowrate' num2str(nc) '.dat'];
    fileID = fopen(filename,'r');
    count = 1;
    values = zeros(nt,3);
    
    while ~feof(fileID)
        values(count,:) = fscanf(fileID, '%f %f %f\n', [1,3]);
        count = count + 1;
    end
    fclose(fileID);
    
    figure
    hold on
    plot(values(1:end-1,1), values(1:end-1,2)/values(end-1,3), 'b-', 'LineWidth', 2)    % plot fracture volume (normalized)
    plot(values(1:end-1,1), values(1:end-1,3)/values(end-1,3), 'r-', 'LineWidth', 2)    % plot total injected volume (normalized)
    hold off
    
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('Normalized Volume [m^3]','FontSize', 20, 'FontWeight', 'bold')
    xlabel('time [s]','FontSize', 20, 'FontWeight', 'bold')
    legend('V_{HF}', 'V_{inj}')
    ylim([0,max(values(:,2))/values(end,3)])
  
    figure
    hold on
    plot(values(1:end,1), values(1:end,2)./values(1:end,3)-1, 'b-', 'LineWidth', 2)
    hold off
    
    title(['Fracture' num2str(nc)], 'FontSize', 20)
    ylabel('error in injected volume','FontSize',20,'FontWeight', 'bold')
    xlabel('time (s)','FontSize',20,'FontWeight', 'bold')
    
end
