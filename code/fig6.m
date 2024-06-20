function fig6

experiment = 1;
close all;

nr = 3;
nc = 2;
fsiz = [0 0 .43 1];
subplots = 1:6;  

figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
hpl_assignment(experiment, nr, nc, subplots(1:4));

model_neutral_response_time(experiment, nr, nc, subplots(5:6));
h(1) = gcf;

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
