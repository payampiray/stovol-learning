function stats = kf_maladaptivity(experiment, nr, nc, subplots)
if nargin<1, experiment = 1; end

[~, lr] = model_neutral(experiment);
lr_eff = lr*[-1 -1 1 1; -1 1 -1 1;1 -1 -1 1]';

% fx = get_kalman(data_set);
[fx, labels, kp] = get_kalman(experiment);

for i=1:3 
    [~,p,~,st]=ttest(fx(:,i));
    st.p = p;
%     st.bf10 = bf.ttest(fx(:,i));
%     st.bf01 = 1/st.bf10;
    st.num_negative = mean(fx(:,i)<0);
    stats.(labels{i}).all = st;
    groups = {fx(:,i)>=0, fx(:,i)<0};
    g_label = {'pos', 'neg'};
    for j=1:2
        g = groups{j};
        [~,p,~,st]=ttest(lr_eff(g,:));
        cc =  corr(lr_eff(g,:));
        st.p = p;
        st.labels = {'PE x Sto', 'PE x Vol', 'PE x Sto x Vol'};
        stats.(labels{i}).(g_label{j}) = st;
        Effect1 = meanEffectSize(lr_eff(g,1), Effect="robustcohen");
        Effect2 = meanEffectSize(lr_eff(g,2), Effect="robustcohen");
    end

end

[cr, cp] = corr(fx(:, 1), fx(:,2), 'type', 'Spearman');
stats.corr_stovol.r = cr;
stats.corr_stovol.p = cp;

neg_eff = sum(fx<0);
q_eff = quantile(fx, [0.25, 0.5, 0.75]);
tbl_data = [q_eff; neg_eff];
stats.table.rows = {'25% quantile'; 'Median'; '75% quantile'; 'Negative value %'};
stats.table.columns = labels;
stats.table.data = tbl_data;


%--------------------------------------------------------------------------
mk = median(kp);
ek = se_median(kp);
mk = [mk([1 3])' mk([2 4])']';
ek = [ek([1 3])' ek([2 4])']';

ma = cell(2,2);
ea = cell(2,2);

s = fx(:, 1)>0;
ma{1,1} = mean(lr(s,:));
ea{1,1} = serr(lr(s,:));
ma{1,2} = mean(lr(~s,:));
ea{1,2} = serr(lr(~s,:));

mf{1}(1,:) = mean(lr(s,:)*[-1 -1 1 1; -1 1 -1 1]');
ef{1}(1,:) = serr(lr(s,:)*[-1 -1 1 1; -1 1 -1 1]');
mf{1}(2,:) = mean(lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]');
ef{1}(2,:) = serr(lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]');
xf{1}{1} = lr(s,:)*[-1 -1 1 1; -1 1 -1 1]';
xf{1}{2} = lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]';

x{1} = {lr(s, [1 3]), lr(s, [2 4])};
x{2} = {lr(~s, [1 3]), lr(~s, [2 4])};
group_title = cell(1,4);
group_title(1:2) = {'$\lambda_s>0$', '$\lambda_s<0$'};


s = fx(:, 2)>0;
ma{2,1} = mean(lr(s,:));
ea{2,1} = serr(lr(s,:));
ma{2,2} = mean(lr(~s,:));
ea{2,2} = serr(lr(~s,:));
group_title(3:4) = {'$\lambda_v>0$', '$\lambda_v<0$'};
x{3} = {lr(s, [1 3]), lr(s, [2 4])};
x{4} = {lr(~s, [1 3]), lr(~s, [2 4])};

mf{2}(1,:) = mean(lr(s,:)*[-1 -1 1 1; -1 1 -1 1]');
ef{2}(1,:) = serr(lr(s,:)*[-1 -1 1 1; -1 1 -1 1]');
mf{2}(2,:) = mean(lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]');
ef{2}(2,:) = serr(lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]');
xf{2}{1} = lr(s,:)*[-1 -1 1 1; -1 1 -1 1]';
xf{2}{2} = lr(~s,:)*[-1 -1 1 1; -1 1 -1 1]';

mx = cell(size(ma));
ex = cell(size(ea));
% x = cell(size(ea));
for i =1:numel(mx)
    mx{i} = [ma{i}([1 3])' ma{i}([2 4])']';
    ex{i} = [ea{i}([1 3])' ea{i}([2 4])']';    
end

%--------------------------------------------------------------------------
if nargin< 2
    close all;    
    nr = 2;
    nc = 3;
    fsiz = [0 0 .65 .55];
    subplots = 1:6;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

colstrs = {'Small','Large'};
glbl = {'Small','Large'};
yl = [0.5 0.82];
fn = def('fn');
fsy = def('fsy');
fst = fsy+4;
yt = 1.05;
lgtitle = def('sto');

elbl = {def('sto'),def('vol')};
egtitle = 'Effect';

% colmap = [255 217 50;97 216 54]/255;
col_br = def('col_yg');

h(1:2) = plot_bar(nr,nc,subplots(1:2), mx(1,:), ex(1,:),colstrs,repmat({'Learning rate coefficient'},1,2));
h(3:4) = plot_bar(nr,nc,subplots(4:5), mx(2,:), ex(2,:),colstrs,repmat({'Learning rate coefficient'},1,2));
lg = legend(h(4),glbl,'fontsize',fsy,'location','north','box','off','AutoUpdate','off');
title(lg,lgtitle,'fontweight','normal');
for i=1:4
    if i>2
        xlabel(h(i),def('vol'),'fontsize',fsy); 
    end
%     text(0.5, yt, group_title{i},'Interpreter','latex','HorizontalAlignment','center','fontsize',fst,'Unit','normalized','fontname',fn,'parent',h(i));
    title(h(i),group_title{i},'Interpreter','latex','fontsize',fst);
    set(h(i),'ylim',yl);
%     plot_dots(h(i), x{i});
end

hh = nan(1,6);
hh([1 2 4 5]) = h;
hh(3) = plot_bar(nr,nc,subplots(3), mf(1), ef(1),group_title(1:2),{'Effect size'},[],col_br);
hh(6) = plot_bar(nr,nc,subplots(6), mf(2), ef(2),group_title(3:4),{'Effect size'},[],col_br);

ii = [3 6];
for i=1:2
    xaxes = get(hh(ii(i)),'XAxis');
    set(xaxes,'TickLabelInterpreter','latex','fontsize',fst);
    plot_dots(hh(ii(i)), xf{i}, false, col_br);
    set(hh(ii(i)),'ylim',[-1.2 1.2]);
end

lg = legend(hh(3),elbl,'fontsize',fsy,'location','south','box','off','AutoUpdate','off','location','southeast');
title(lg,egtitle,'fontweight','normal');


% % -----------------
% % Figure 2
% figure;
% plot_bar(1,1,1, {mk}, {ek}, {'Small','Large'},{'Kalman parameter'});
end

function [fx, labels, kalman_parameter] = get_kalman(data_set)
fname = fullfile('..',sprintf('experiment%d', data_set),'model_kf4.mat');

if ~exist(fname, 'file')
    data = get_data(data_set);    
    num_params = 4;
    N = length(data);
    lme = nan(N, 1);
    effect = nan(N, num_params);
    kalman_parameter = nan(N, num_params);
    dynamics = cell(N, 1);
    for n=1:N        
        kalman(n) = fit1_fit(data{n}, []); %#ok<AGROW> 
        lme(n) = kalman(n).lme;
        effect(n, :) = kalman(n).effect;
        kalman_parameter(n, :) = kalman(n).kalman_parameter;
        dynamics{n} = kalman_model(kalman(n).kalman_parameter, data{n}.bag);
        fprintf('%03d\n', n);
    end

    save(fname, 'kalman', 'lme', 'effect', 'kalman_parameter', 'dynamics');
end

f = load(fname);
fx = f.effect;
fx = fx(:, [2:4 1]);
labels = {'lambda_s', 'lambda_v', 'lambda_i', 'lambda_m'};
kalman_parameter = f.kalman_parameter;
end

function [dynamics] = kalman_model(lambda, outcome)
ndim = size(outcome,2);
N = size(outcome,1);
lr = nan(N,ndim);
val = nan(N,ndim);

m = 60+zeros(1,ndim);
r = 1*ones(1,ndim);

for t=1:N
    val(t, :) = m;
    a = (r+lambda)./(r+1+lambda);
    m = m + a.*(outcome(t, :)-m);
    r = (1-a).*(r+lambda); 

    lr(t, :) = a;
end

dynamics = struct('val', val, 'lr', lr);

end

