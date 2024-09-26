function fig6

experiment = 1;
close all;

[~, h1] = hpl_assignment(experiment, 2, 1, 1:2);

[~, h2] = model_neutral_response_time(experiment, 1, 1, 1);
h = [h1, h2];

for i = 1:length(h)
    figname = fullfile('..','figs', sprintf('%s_%d', mfilename, i));
    saveas(h(i), figname, 'jpg');
end

end
