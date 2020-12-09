
global IOPath

filename = [IOPath 'WB_NR_FC.dat'];
fileID = fopen(filename,'r');
count = 1;

while ~feof(fileID)
    temp = fscanf(fileID, '%f\n', [1,1]);
    if temp ~= 0
        values(count,:) = temp;
        count = count + 1;
    end
end
fclose(fileID);

figure
semilogy(1:length(values), values)
ylabel('L2 norm of error','FontSize', 24)
xlabel('iterations','FontSize', 24)

