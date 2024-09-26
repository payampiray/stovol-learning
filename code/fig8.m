function fig8

experiment = 2;

close all;
nr = 2;
nc = 3;
fsiz = [0 0 .57 .55];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

hpl_clustering(experiment, nr, nc, 1:3, .5);
h(1) = gcf;

nr = 1;
fsiz = [0 0 .57 .225];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_assignment(experiment, nr, nc, 1:2, 2);
model_neutral_response_time(experiment, nr, nc, 3, 2)
h(2) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
