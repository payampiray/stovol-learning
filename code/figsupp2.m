function figsupp2

experiment = 1;

close all;


nr = 2;
nc = 3;
subplots = [1 2];
fsiz = [0 0 .7 .55];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
[~, h1] = hpl_dynamics(experiment, nr, nc, subplots);

% nr = 1;
% nc = 3;
subplots = 3:6;
% fsiz = [0 0 .7 .27];
% figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

[~, h2] = hpl_modulation(experiment, nr, nc, subplots);

h = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
