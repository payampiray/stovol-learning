function fig4

data_set = 1;

close all;

fsiz = [0 0 .85 .3];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_signal(data_set, 1, 4, 2:4, 1);
h(1) = gcf;

fsiz = [0 0 .85 .55];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_signal(data_set, 2, 2, 1:4, 2);
h(2) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end