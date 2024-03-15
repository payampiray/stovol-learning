function stats = hpl_clustering(experiment, nr, nc, subplots, fig_no)
% make_signal
if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end

data = get_data(experiment);

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fsig,'file')
    for n=1:length(data)
        dat = data{n};                
        [a_data(n, :)] = tools_clustering(dat.bucket, dat.bag);

        dynamics = f.dynamics{n};
        [a_model(n, :), s_model(n, :), v_model(n, :), labels] = tools_clustering(dynamics.val, dat.bag, dynamics.sto, dynamics.vol);
    end
    save(fsig, 'a_data', 'a_model', 's_model', 'v_model', 'labels');
end
f = load(fsig);
a_data = f.a_data;
a_model = f.a_model;
s_model = f.s_model;
v_model = f.v_model;
labels = f.labels;

% ------------------
x = {v_model, s_model, a_model, a_data};
for i=1:length(x)
    mdb(i, :) = mean(x{i});
    sdb(i, :) = serr(x{i});
    [~, pdb(i, :),~, tstat(i)] = ttest(x{i});    

    xdb(i, :) = median(x{i});
    edb(i, :) = se_median(x{i}, 5000);
    [qdb(i, :),~, wilcox(i)] = signrank2(x{i});
end

y_labels = {'v_model', 's_model', 'a_model', 'a_data'};
for i=1:4
    st = tstat(i);
    st.p = pdb(i, :);
    st.mean = mdb(i, :);
    st.serr = sdb(i, :);

    st.q = qdb(i, :);
    st.zval = wilcox(i).zval;
    st.median = xdb(i, :);
    st.se_median = edb(i, :);
    
    st.columns = labels;
    stats.(y_labels{i}) = st;
end

%--------------------------------------------------------------------------
mx = mdb(3:4, :);
ex = sdb(3:4, :);
bs = a_data(:, 1);
bv = a_data(:, 2);

if fig_no == .5
    if nargin< 2
        close all;
    
        nr = 2;
        nc = 3;
        fsiz = [0 0 .7 .6];  
        subplots = 1:4;
        
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end
    
    colstrs = labels;
    fs = def('fs');
    fsy = def('fs');
    fst = def('fsy');
    y_str1 = sprintf('Changes in learning rate');
    
    colmap = def('col_yg');
    bw1 = 0.05;
    
    plot_bar(nr, nc, subplots(1), {mx(1,:)}, {ex(1,:)}, colstrs, {y_str1}, [], colmap, {''}, bw1);
    title('Model');
    
    % ar = get(gca,'PlotBoxAspectRatio');
    % ar(2) = .9*ar(2);
    % set(gca, 'PlotBoxAspectRatio', ar);  
    xlabel('Sample Autocorrelation', 'fontsize', fst)
    
    plot_bar(nr, nc, subplots(2), {mx(2,:)}, {ex(2,:)}, colstrs, {y_str1}, [], colmap, {''}, bw1);
    ytick = get(gca, 'ytick');
    yt = cellstr(num2str(ytick'));
    set(gca, 'yticklabel', yt);
    title('Data');
    xlabel('Sample Autocorrelation', 'fontsize', fst);


    return;
end
%--------------------------------------------------------------------------


mx = [xdb(1:2, :); mdb(3:4, :)];
ex = [edb(1:2, :); sdb(3:4, :)];
bs = a_data(:, 1);
bv = a_data(:, 2);


if nargin< 2
    close all;
%     nr = 2;
%     nc = 2;
%     fsiz = [0 0 .55 .65];
%     subplots = 1:4;    

    nr = 2;
    nc = 3;
    fsiz = [0 0 .7 .6];  
    subplots = 1:6;

    nr = 3;
    nc = 2;
    fsiz = [0 0 .4 .9];  
    subplots = 1:6;    
    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

colstrs = labels;
fs = def('fs');
fsy = def('fs');
fst = def('fsy');
y_str1 = sprintf('Changes in learning rate');
y_str2 = sprintf('Changes in learning rate');
y_str = {sprintf('Changes in \nvolatility estimate'), sprintf('Changes in \nstochasticity estimate'), ...
         'Changes in learning rate', 'Changes in learning rate'};
titles = {'Model', 'Model', 'Model', 'Data'};

colmap = def('col_yg');
bw1 = 0.05;

for i=1:4
    plot_bar(nr, nc, subplots(i), {mx(i,:)}, {ex(i,:)}, colstrs, y_str(i), [], colmap, {''}, bw1);
    title(titles{i});
    xlabel('Sample Autocorrelation', 'fontsize', fst)
end

% ---------------------------------
% 2nd figure

subplot(nr, nc, subplots(5:6));

rand_sigma = 0.05;
rng(0);
for n=1:size(bs,1)
    plot([bs(n), bv(n)], [1 2]+randn*rand_sigma, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
end

config.dens_ratio = 2.5;
config.face_color = colmap;
config.patch_alpha = .4;
y_points = [0.65 2.35];
raincloud1xN_horizontal([bs bv], y_points, config)

h = set(gca, 'ytick', y_points, 'yticklabel', colstrs, 'fontsize', fsy);
title('Data');
xlabel(y_str2, 'fontsize', fst);
ax = ancestor(gca, 'axes');
xaxes = get(ax,'YAxis');
set(xaxes, 'TickLabelRotation', 90, 'fontsize', fst);
ylabel('Sample Autocorrelation', 'fontsize', fst);
yt = get(gca,'ylim');
set(gca,'ylim', [-.7 3.5])
end
