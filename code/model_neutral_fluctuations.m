function model_neutral_fluctuations(experiment, nr, nc, subplots)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('model_neutral_fluctuations.mat'));
if 1 %~exist(fname, 'file')
    [data] = get_data(experiment);
    [x, e] = tools_fluctuations(data);
    save(fname, 'x', 'e');
end
f = load(fname);
x = f.x;
e = f.e;

if nargin< 2
    close all;    
    nr = 1;
    nc = 2;
    fsiz = [0 0 .65 .25];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

fsy = def('fsy');
fs = def('fs');
cols = def('col');

ii = {[1 3], [2 4]};
lg_title = def('sto');
fig_titles = {sprintf('Small %s', lower(def('vol'))), sprintf('Large %s', lower(def('vol')))};

% ii = {[1 2], [3 4]};
% lg_title = def('vol');
% fig_titles = {sprintf('Small %s', lower(def('sto'))), sprintf('Large %s', lower(def('sto')))};

label = 'Learning rate coefficient';
lg_labels = {'Small', 'Large'};
yl = [.2 1.3];
for i=1:2
    subplot(nr, nc, subplots(i));
    for j=1:2
    hm(j) = plot_ci(x(:, ii{i}(j)), e(:, ii{i}(j)), cols(j, :)); hold on;
    ylim(yl);
    end
    ylabel(label, 'fontsize', fsy);
    xlabel('Trial', 'fontsize', fsy);
    title(fig_titles{i}, 'fontsize', fsy, 'fontweight', 'normal');
    set(gca, 'box', 'off');

    if i==2
        lg = legend(hm, lg_labels, 'box', 'off', 'location', 'southeast', 'fontsize', fs);
        title(lg, lg_title,'fontweight','normal', 'fontsize', fsy);
    end
end

end
