function fig8

data_set = 2;

close all;
nr = 3;
nc = 2;
fsiz = [0 0 .4 .85];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

hpl_signal(data_set, nr, nc, 1:2, .5);

hpl_clustering(data_set, nr, nc, 3:4, .5);

hpl_modulation(data_set, nr, nc, 5:6, .5);


h(1) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
