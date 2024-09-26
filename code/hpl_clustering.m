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
    [~, pdb(i, :), ci{i}, tstat(i)] = ttest(x{i});    

    xdb(i, :) = median(x{i});
    edb(i, :) = se_median(x{i}, 5000);
    [qdb(i, :),~, wilcox(i)] = signrank2(x{i});
end

y_labels = {'v_model', 's_model', 'a_model', 'a_data'};
for i=1:4
    st = tstat(i);
    st.p = pdb(i, :);
    st.ci = ci{i};
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

%--------------------------------------------------------------------------


mx = [xdb(1:2, :); mdb(3:4, :)];
ex = [edb(1:2, :); sdb(3:4, :)];

bs(:, 1) = a_model(:, 1);
bv(:, 1) = a_model(:, 2);

bs(:, 2) = a_data(:, 1);
bv(:, 2) = a_data(:, 2);


colstrs = labels;
fs = def('fs');
fsy = def('fs')+2;
fst = def('fsy');
y_str1 = sprintf('Changes in learning rate');
y_str2 = sprintf('Changes in learning rate');
y_str = {sprintf('Changes in \nvolatility estimate'), sprintf('Changes in \nstochasticity estimate'), ...
         'Changes in learning rate', 'Changes in learning rate'};
titles = {'Model', 'Data'};
colmap = def('col_yg');
bw1 = 0.05;

if fig_no == .5
    if nargin< 2
        close all;
    
        nr = 1;
        nc = 3;
        fsiz = [0 0 .5 .25];  
        subplots = [1 2 3];    
        
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end

    mx = mx(4, :);
    ex = ex(4, :);
    y_str = y_str(4);
    bs = bs(:, 2);
    bv = bv(:, 2);

    for i=1
        plot_bar(nr, nc, subplots(1), {mx(i,:)}, {ex(i,:)}, colstrs, y_str(i), [], colmap, {''}, bw1);
        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'XAxis');
        set(xaxes,'fontsize',fsy);  
        xlabel('Sample Autocorrelation', 'fontsize', fst);        
    end
    
    for i=1
        subplot(nr, nc, subplots(2:3));        
        rand_sigma = 0.05;
        rng(0);
        for n=1:size(bs,1)
            plot([bs(n, i), bv(n, i)], [1 2]+randn*rand_sigma, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
        end
        
        config.dens_ratio = 2.5;
        config.face_color = colmap;
        config.patch_alpha = .4;
        y_points = [0.65 2.35];
        raincloud1xN_horizontal([bs(:, i) bv(:, i)], y_points, config)
        
        h = set(gca, 'ytick', y_points, 'yticklabel', colstrs, 'fontsize', fs);
        xlabel(y_str2, 'fontsize', fst);
        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'YAxis');
        set(xaxes, 'TickLabelRotation', 90, 'fontsize', fsy);
        ylabel('Sample Autocorrelation', 'fontsize', fst);
        yt = get(gca,'ylim');
        set(gca,'ylim', [-1 4])

    end

    return;
end
%--------------------------------------------------------------------------


mx = [xdb(1:2, :); mdb(3:4, :)];
ex = [edb(1:2, :); sdb(3:4, :)];

bs(:, 1) = a_model(:, 1);
bv(:, 1) = a_model(:, 2);

bs(:, 2) = a_data(:, 1);
bv(:, 2) = a_data(:, 2);


colstrs = labels;
fs = def('fs');
fsy = def('fs')+2;
fst = def('fsy');
y_str1 = sprintf('Changes in learning rate');
y_str2 = sprintf('Changes in learning rate');
y_str = {sprintf('Changes in \nvolatility estimate'), sprintf('Changes in \nstochasticity estimate'), ...
         'Changes in learning rate', 'Changes in learning rate'};
titles = {'Model', 'Data'};

colmap = def('col_yg');
bw1 = 0.05;

N = size(v_model, 1);
if fig_no == 1
    if nargin< 2
        close all;
    
        nr = 1;
        nc = 2;
        fsiz = [0 0 .5 .2];  
        subplots = [1 2];    
        
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end

    for i=1:2
    
        h(i) = subplot(nr, nc, subplots(i));
        index = [];
        data = [];
        for j=1:size(x{i}, 2)
            data = [data; x{i}(:, j)];
            index = [index; j*ones(size(x{i}, 1), 1)];
        end
        tbl = table(data);
        tbl.index = index;
        box = boxchart(tbl.index, tbl.data, 'GroupByColor',tbl.index);
        for j=1:length(box)
            set(box(j), 'BoxEdgeColor', [0 0 0], 'BoxFaceColor', colmap(j, :), 'BoxMedianLineColor', [colmap(j, :) ], 'MarkerColor', colmap(j, :), 'BoxFaceAlpha', .6);        
    
        end
        set(gca, 'xtick', [.75 2.25], 'xticklabel',colstrs, 'xlim', [0 3], 'fontsize', fs);
        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'XAxis');
        set(xaxes,'fontsize',fsy);      
    
    
        set(gca,'box','off');
        alpha(.6);    
    % 
        title('Model', 'fontsize', fst);
        xlabel('Sample Autocorrelation', 'fontsize', fst);
        ylabel(y_str{i}, 'fontsize', fst);
    end
end

% ---------------------------------
if fig_no == 2
    if nargin< 2
        close all;
    
        fsiz = [0 0 .55 .55];          
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end
    
    nr = 2;
    nc = 3;    
    subplots = {1, 4, 2:3, 5:6};    

    mx = mx(3:4, :);
    ex = ex(3:4, :);
    y_str = y_str(3:4);
    for i=1:2
        plot_bar(nr, nc, subplots{i}, {mx(i,:)}, {ex(i,:)}, colstrs, y_str(i), [], colmap, {''}, bw1);
        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'XAxis');
        set(xaxes,'fontsize',fsy);  
        xlabel('Sample Autocorrelation', 'fontsize', fst);
        title(titles{i}, 'fontsize', fst);
    end
    
    for i=1:2
        subplot(nr, nc, subplots{i+2});        
        rand_sigma = 0.05;
        rng(0);
        for n=1:size(bs,1)
            plot([bs(n, i), bv(n, i)], [1 2]+randn*rand_sigma, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
        end
        
        config.dens_ratio = 2.5;
        config.face_color = colmap;
        config.patch_alpha = .4;
        y_points = [0.65 2.35];
        raincloud1xN_horizontal([bs(:, i) bv(:, i)], y_points, config)
        
        h = set(gca, 'ytick', y_points, 'yticklabel', colstrs, 'fontsize', fs);
        title(titles{i}, 'fontsize', fst);
        xlabel(y_str2, 'fontsize', fst);
        ax = ancestor(gca, 'axes');
        xaxes = get(ax,'YAxis');
        set(xaxes, 'TickLabelRotation', 90, 'fontsize', fsy);
        ylabel('Sample Autocorrelation', 'fontsize', fst);
        yt = get(gca,'ylim');
        set(gca,'ylim', [-1 4])

    end
end
