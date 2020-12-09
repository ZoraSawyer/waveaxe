global IOPath

filename = [IOPath 'WB_NR_UD.dat'];
fileID = fopen(filename,'r');
count = 1;
niter = 200;
values = zeros(niter,1);
plot_freq = [1,16:20:96];
hn = 1:length(plot_freq);

while ~feof(fileID)
    values(:,count) = fscanf(fileID, '%e\n', [niter,1]);
    count = count + 1;
end
fclose(fileID);

for n = 1:length(plot_freq)
    index = find(values(:,plot_freq(n)) ~= 0);
    v = values(index,plot_freq(n));
    figure(hn(n))
    hold on
    semilogy(1:length(v), v)
    hold off
    ylabel('error','FontSize', 24)
    xlabel('iteration','FontSize', 24)
end
