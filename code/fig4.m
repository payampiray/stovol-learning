function fig4

experiment = 1;

close all;

fsiz = [0 0 .18 .22];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_signal(experiment, 1, 1, 1, 1);
h(1) = gcf;

fsiz = [0 0 .3 .45];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_signal(experiment, 2, 2, 1:2, 1.5);
h(2) = gcf;

fsiz = [0 0 .7 .45];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_signal(experiment, 2, 2, 1:4, 2);
h(3) = gcf;


end