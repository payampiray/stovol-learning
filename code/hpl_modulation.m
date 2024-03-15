function stats = hpl_modulation(experiment, nr, nc, subplots, fig_no)
% make_signal
if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end

data = get_data(experiment);

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fsig,'file')
    for n=1:length(data)
        dat = data{n};                
        [b_data(n, :), b_labels] = tools_modulation(dat.bucket, dat.bag);

        dynamics = f.dynamics{n};
        [b_model(n, :)] = tools_modulation(dynamics.val, dat.bag);
    end
    save(fsig, 'b_data', 'b_labels', 'b_model');
end
f = load(fsig);
b_data = f.b_data;
b_model = f.b_model;
b_labels = f.b_labels;

x = {b_model, b_data};
for i=1:length(x)
    x{i} = x{i};
    xdb(i, :) = median(x{i});
    mdb(i, :) = mean(x{i});
    sdb(i, :) = serr(x{i});
    edb(i, :) = se_median(x{i});
    [qdb(i, :),~, wilcox(i)] = signrank2(x{i});
    [~, pdb(i, :),~, tstat(i)] = ttest(x{i});
    tval(i, :) = tstat(i).tstat;
    se(i, :) = tstat(i).sd;
    zval(i, :) = wilcox(i).zval;
end

% lr = glm1(data_set);
% fx = lr*[-1 -1 1 1;-1 1 -1 1]';
% fx = fx.*[-1 1];
% fx(:, 3) = 2*(sum(fx, 2)>0)-1;
% dlr = x{2};
% for i=1:3
%     mdlr_eff(i, :) = mean(dlr(fx(:, i)>0, :));
%     [~, s(i,:)] = ttest2(dlr(fx(:,i)>0, :), dlr(fx(:,i)<0, :));
%     [sx(i,:)] = ranksum2(dlr(fx(:,i)>0, :), dlr(fx(:,i)<0, :));
%     md(i,:) = mean(dlr(fx(:,i)<0, :)) - mean(dlr(fx(:,i)>0, :));
%     xd(i,:) = median(dlr(fx(:,i)<0, :)) - median(dlr(fx(:,i)>0, :));
% end

y_labels = {'model', 'data'};
i = 1;
st = wilcox(i);
st.p = qdb(i, :);
st.median = xdb(i, :);
st.se_median = edb(i, :);
st.columns = b_labels;
stats.(y_labels{i}) = st;

i = 2;
st = tstat(i);
st.p = pdb(i, :);
st.mean = mdb(i, :);
st.se = sdb(i, :);
st.columns = b_labels;
stats.(y_labels{i}) = st;

tbl_data = [xdb(1, :); zval(1,:); qdb(1, :); mdb(2, :); sdb(2, :); tval(2, :); pdb(2, :)];

% stats.table(1).rows = {'Median for model LR'; 'z-val for model LR'; 'P-value for model LR'};
% stats.table(1).columns = b_labels;
% stats.table(1).data = tbl_data(1:3, :);

stats.table(1).columns = {'Estimate for model-agnostic LR'; 'SE for model-agnostic LR'; 't-value for model LR'; 'P-value for model LR'};
stats.table(1).rows = b_labels;
stats.table(1).data = tbl_data(4:end, :)';


%--------------------------------------------------------------------------

mx = mdb(:, 5:6);
ex = sdb(:, 5:6);
bs = b_data(:, 5);
bv = b_data(:, 6);


if nargin< 2
    close all;

    nr = 2;
    nc = 3;
    fsiz = [0 0 .7 .6];
    subplots = 3:6;    
    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

% colstrs = {'$|\delta|$', '$|\Delta m|$'};
colstrs = {'|PE|+|PU|', '|PE|-|PU|'};
% colstrs = {'Small','Large'};
fn = def('fn');
fs = def('fs');
fsy = def('fsy');
fst = fsy;
y_str1 = sprintf('Modulation of learning rate\n (regression coefficient)');
y_str2 = sprintf('Modulation of learning rate (regression coefficient)');
y_str1 = sprintf('Modulation of learning rate');
y_str2 = sprintf('Modulation of learning rate');


colmap = def('col_yg');
bw1 = 0.05;

plot_bar(nr,nc, subplots(1), {mx(1,:)}, {ex(1,:)}, colstrs, {y_str1}, [], colmap, {''}, bw1);
title('Model');

ar = get(gca,'PlotBoxAspectRatio');
ar(2) = .9*ar(2);
set(gca, 'PlotBoxAspectRatio', ar);    

plot_bar(nr,nc,subplots(2), {mx(2,:)}, {ex(2,:)}, colstrs, {y_str1}, [], colmap, {''}, bw1);
ytick = get(gca, 'ytick');
yt = cellstr(num2str(ytick'));
set(gca, 'yticklabel', yt);
title('Data');

ar = get(gca,'PlotBoxAspectRatio');
ar(2) = .9*ar(2);
set(gca, 'PlotBoxAspectRatio', ar);

if fig_no == 0.5
    return; 
end

% ---------------------------------
% 2nd figure

subplot(nr, nc, subplots(3:4));

rand_sigma = 0.05;
rng(0);
for n=1:size(bs,1)
    plot([bv(n), bs(n)], [1 2]+randn*rand_sigma, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
end

colmap = colmap([2 1], :);
colstrs = colstrs([2 1]);

config.dens_ratio = 2.5;
config.face_color = colmap;
config.patch_alpha = .4;
y_points = [0.65 2.35];
raincloud1xN_horizontal([bv bs], y_points, config)

set(gca, 'ytick', y_points, 'yticklabel', colstrs, 'fontsize', fs);
title('Data');
xlabel(y_str2)

% rand_sigma = 0.05;
% rng(0);
% for n=1:size(bs,1)
%     plot([bv(n)], 2+randn*rand_sigma, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
% end
% 
% config.dens_ratio = 2.5;
% config.face_color = colmap(2,:);
% config.patch_alpha = .4;
% y_points = [2.35];
% raincloud1xN([bv], y_points, config)
% 
% set(gca, 'ytick', y_points, 'yticklabel', colstrs, 'fontsize', fs);
% title('Model-agnostic learning rate');
% xlabel(y_str)

end

