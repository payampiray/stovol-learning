function fig3

data_set = 1;

nr = 2;
nc = 3;
fsiz = [0 0 .65 .55];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

kf_maladaptivity(data_set, nr,nc, 1:6);
h(1) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end


end
