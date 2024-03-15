function figsupp1


nr = 1;
nc = 3;
subplots = 1:3;
fsiz = [0.1    0.0800    .75    .3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

hpl_recovery(nr, nc, subplots);


end
