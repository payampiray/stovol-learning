function fig5

data_set = 1;

nr = 1;
nc = 2;
fsiz = [0 0 .6 .2];

close all

figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_clustering(data_set, nr, nc, 1:2, 1);
h(1) = gcf;

fsiz = [0 0 .6 .5];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_clustering(data_set, 2, 2, 1:2, 2);
h(2) = gcf;

end

