function fig1


fsiz = [0 0 .3 .45];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 2;
nc = 1;
design(nr,nc,1:2, 1);
h(1) = gcf;


fsiz = [0 0 .4 .25];
nr = 1;
nc = 2;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
design(nr,nc,1:2, 2);
h(2) = gcf;

fsiz = [0 0 .18 .25];
nr = 1;
nc = 1;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
design(nr, nc, 1, 3);
h(3) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
