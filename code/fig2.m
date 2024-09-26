function fig2

data_set = 1;

fsiz = [0 0 .18 .3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

model_neutral(data_set, 1,1, 1, 1);
h(1) = gcf;

fsiz = [0 0 .4 .3];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
model_neutral(data_set, 1,2, 1:2, 2);
h(2) = gcf;

fsiz = [0 0 .65 .25];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
model_neutral_fluctuations(data_set, 1, 2, 1:2);
h(3) = gcf;


end
