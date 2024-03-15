function fig5

data_set = 1;

nr = 3;
nc = 2;
fsiz = [0 0 .4 1];  
subplots = 1:6;  

close all

figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_clustering(data_set, nr, nc, subplots);
h(1) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
