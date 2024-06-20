function fig8

experiment = 2;

close all;
nr = 1;
nc = 2;
fsiz = [0 0 .4 .2];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

hpl_signal(experiment, nr, nc, 1:2, .5);
h(1) = gcf;

figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_clustering(experiment, nr, nc, 1:2, .5);
h(2) = gcf;

nr = 2;
fsiz = [0 0 .4 .4];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_assignment(experiment, nr, nc, 1:2, 0.5);
model_neutral_response_time(experiment, nr, nc, 3:4)
h(3) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
