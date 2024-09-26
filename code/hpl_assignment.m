function [stats, figs] = hpl_assignment(experiment, nr, nc, subplots, fig_no)

if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end
[data] = get_data(experiment);

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fsig, 'file')
    for n=1:length(data)
        dat = data{n};                
        [b_data(n, :), x_data(n,:)] = tools_assignment(dat.bucket, dat.bag);
    
        dynamics = f.dynamics{n};
        [b_model(n, :), x_model(n, :), bv2s(n, :)] = tools_assignment(dynamics.val, dat.bag, dynamics.vol, dynamics.sto);
    
        save(fsig, 'b_data', 'x_data', 'b_model', 'x_model', 'bv2s');
    end
end

f = load(fsig);
b_data = f.b_data;
x_data = f.x_data;
b_model = f.b_model;
x_model = f.x_model;
bv2s = f.bv2s;

b = {b_model, b_data, bv2s};
for i=1:length(b)
    xdb(i, :) = median(b{i});
    mdb(i, :) = mean(b{i});
    sdb(i, :) = serr(b{i});
    [qdb(i, :),~, wilcox(i)] = signrank2(b{i});
    [~, pdb(i, :), ci{i}, tstat(i)] = ttest(b{i});
    tval(i, :) = tstat(i).tstat;
    se(i, :) = tstat(i).sd;
end

x = {x_model, x_data};
for i=1:length(x)
    mx(i, :) = mean(x{i});
    sx(i, :) = serr(x{i});
end



y_labels = {'model', 'data', 'b2v'};
for i=1:3
    st = tstat(i);
    st.p = pdb(i, :);
    st.ci = ci{i};
    st.mean = mdb(i, :);
    st.sem = sdb(i, :);
    st.columns = {'|AC|', 'AC', 'intercept'};
    stats.(y_labels{i}) = st;
end

%--------------------------------------------------------------------------
mb = mdb(1:3, 1);
eb = sdb(1:3, 1);

b_dots = [b_model(:, 1), b_data(:, 1), bv2s(:, 1)];
%--------------------------------------------------------------------------
if fig_no == 2
    if nargin< 2
        close all;
        nr = 1;
        nc = 2;
        fsiz = [0 0 .4 .2];  
        subplots = [1 2];
        
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end
    figs(1) = gcf;
    
    fsy = def('fsy');
    yl = [-19 9]*10^-4;
    
    yls{1} = sprintf('Relationship between\n |LR| changes and |AC|');
    yls{2} = sprintf('Relationship between\n |LR| changes and |AC|');
    yls{3} = sprintf('Relationship between\nlog %s and |AC|', '\DeltaV/\DeltaS');
    title_str = {'Model', 'Data', 'Model'};
    
    ii = 1:2;
    
    for i=ii
        subplot(nr, nc, subplots(i));    
        plot_raincloud(mb(i), eb(i), b_dots(:, i), experiment);
        ylabel(yls{i}, 'fontsize', fsy);
        title(title_str{i}, 'fontsize', fsy);
        ylim(yl);
    
    end
    return;
end

%--------------------------------------------------------------------------

if nargin< 2
    close all;
    nr = 2;
    nc = 1;
    subplots = [1, 2];
end

fsy = def('fsy');
fs = def('fs');
colmap = def('col_unique');
bw1 = 0.08;

yls{1} = sprintf('Relationship between\n |LR| changes and |AC|');
yls{2} = sprintf('Relationship between\n |LR| changes and |AC|');
yls{3} = sprintf('Relationship between\nlog %s and |AC|', '\DeltaV/\DeltaS');
title_str = {'Model', 'Data', 'Model'};

xl = 'Relationship between |LR| changes and |AC|';

ii = 1:2;

fsiz = [0 0 .16 .5];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);     
for i=1:2
%     plot_bar(nr, nc, subplots(1), {mx(i,:)}, {eb(i,:)}, {''});
    plot_bar(nr, nc, subplots(i), {mb(i)}, {eb(i)}, {''}, yls(1), [], colmap, {''}, bw1);
    title(title_str{i}, 'fontsize', fsy);        
end    
figs(1) = gcf;

fsiz = [0 0 .3 .5];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);  
for i=ii
    subplot(nr, nc, subplots(i));    
%     plot_raincloud(mb(i), eb(i), b_dots(:, i), experiment);

    config.add_dots = 1;
    config.is_right = 1;
    config.dens_ratio = 2.5;
    config.face_color = colmap;
    config.patch_alpha = .4;
    config.func_summary = 'median';
    y_points = [0.5 3.5];
    raincloud1xN_horizontal(b_dots(:, i), y_points, config)    
    set(gca, 'fontsize', fs, 'ytick', []);

    xlabel(xl, 'fontsize', fsy);
    title(title_str{i}, 'fontsize', fsy);    
%     ylim(yl);

end
figs(2) = gcf;

fsiz = [0 0 .2 .5];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);  

yls = {[0 0.2],[.2 .6]};
cols = repmat(colmap(1, :), 10, 1);

wbin = 100/length(mx);
xperc = wbin:wbin:(100);
y_str = 'Changes in |LR|';
% plot_bar(nr, nc, subplots(4), {mx}, {ex}, '', {y_str}, [], cols, {''}, bw1, 1);

for i=ii
    subplot(nr, nc, subplots(i));
    shadedErrorBar(xperc, mx(i, :), sx(i, :), 'lineProps', {'','color',cols(1,:), 'markerfacecolor', [0 0 0]}); hold on;
    plot(xperc, mx(i, :), '.','markerfacecolor', [0 0 0], 'markersize', 10); hold on;


    set(gca, 'fontsize', fs);
    ylabel(y_str, 'fontsize', fsy);
    xlabel('|AC|% (binned)', 'fontsize', fsy);
    ylim(yls{i});
    title(title_str{i}, 'fontsize', fsy);    
end
figs(3) = gcf;

end
