function fig6

data_set = 1;

close all;


nr = 2;
nc = 3;
subplots = [1 2];
fsiz = [0 0 .7 .55];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_dynamics(data_set, nr, nc, subplots);
% h(1) = gcf;

% nr = 1;
% nc = 3;
subplots = 3:6;
% fsiz = [0 0 .7 .27];
% figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

hpl_modulation(data_set, nr, nc, subplots);
h(1) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
