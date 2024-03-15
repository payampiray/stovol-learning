function fig7

data_set = 2;

nr = 1;
nc = 1;
fsiz = [0 0 .18 .25];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
model_neutral(data_set, nr,nc,1, 1);
h(1) = gcf;

fsiz = [0 0 .4 .25];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
model_neutral(data_set, nr,nc,1, 2);
h(2) = gcf;

nr = 2;
nc = 3;
fsiz = [0 0 .65 .55];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
kf_maladaptivity(data_set, nr,nc, 1:6);
h(3) = gcf;

for i = 1:length(h)
    figname = fullfile('figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end


end
