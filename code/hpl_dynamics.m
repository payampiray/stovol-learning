function [stats, h] = hpl_dynamics(experiment, nr, nc, subplots, fig_no)

if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end

data = get_data(experiment);

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fsig,'file')
    N = length(data);
    for n=1:N
        dat = data{n};
        dynamics = f.dynamics{n};
        [beta_vol(n, :), beta_sto(n, :), beta_difr(n,:), x_labels, y_labels] = regression(dat.bag, dynamics.val, dynamics.vol, dynamics.sto);
    end
    save(fsig, 'beta_vol', 'beta_sto', 'beta_difr', 'x_labels', 'y_labels');
end
f = load(fsig);
beta_vol = f.beta_vol;
beta_sto = f.beta_sto;
beta_difr = f.beta_difr;
x_labels = f.x_labels;
y_labels = f.y_labels;

x = {beta_vol, beta_sto, beta_difr};
for i=1:length(x)
    x{i} = x{i};
    xdb(i, :) = median(x{i});
    mdb(i, :) = mean(x{i});
    sdb(i, :) = serr(x{i});
    edb(i, :) = se_median(x{i});
    [qdb(i, :),~, wilcox(i)] = signrank2(x{i});
    zdb(i,:) = wilcox(i).zval;
    [~, pdb(i, :),~, stt(i)] = ttest(x{i});

    ndb(i, :) = mean(x{i}>0);
end

for i=1:length(y_labels)
    st = wilcox(i);
    st.median = mdb(i, :);
    st.p = qdb(i, :);
    st.columns = x_labels;
    stats.(y_labels{i}) = st;
end
% 
tbl_data = [zdb; qdb];
stats.table.rows = {
                 'z-value for volatility dynamics'; 'z-value for stochasticity dynamics'; 'z-value for the log ratio';...
                 'P-value for volatility dynamics'; 'P-value for stochasticity dynamics'; 'P-value for the log ratio'};
stats.table.columns = x_labels;
stats.table.data = tbl_data;

%--------------------------------------------------------------------------

mx = xdb(:, 2:3);
ex = edb(:, 2:3);


% colstrs = {'$|\delta|$', '$|\Delta m|$'};
colstrs = {'|PE|+|PU|', '|PE|-|PU|'};
% colstrs = {'Small','Large'};
fn = def('fn');
fs = def('fs');
fsy = def('fsy');
fst = fsy;
title_str = {sprintf('Trial-by-trial dynamics'),sprintf('Volatility/Stochasticity'),'Model learning rate'};
% ylabel_str = {sprintf('Trial-by-trial dynamics\nregression coefficient'),...
%               sprintf('Dynamics of V/S ratio\n regression coefficient'),...
%               sprintf('Modulation of learning rate\n regression coefficient')};
ylabel_str = {sprintf('Trial-by-trial dynamics'),...
              sprintf('Dynamics of V/S ratio'),...
              sprintf('Modulation of learning rate')};
title_str = repmat(title_str, 1, 2);
y_str = 'Regression coefficient';

colmap = def('col_yg');
bw1 = 0.1;
bw2 = 0.05;

if fig_no == 1
    if nargin< 2
        close all;    
        nr = 1;
        nc = 3;
        fsiz = [0 0 .7 .3];
        subplots = 1:3;
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end
    
    h(1) = plot_bar(nr,nc,subplots(1), {mx(1:2,:)}, {ex(1:2,:)},{'Volatility', 'Stochasticity'}, {y_str}, [], colmap, {''}, bw1);
    h(2) = plot_bar(nr,nc,subplots(2), {mx(3,:)}, {ex(3,:)}, colstrs, {y_str}, [], colmap, {''}, bw2);
%     h(3) = plot_bar(nr,nc,subplots(3), {mx(4,:)}, {ex(4,:)}, colstrs, {y_str}, [], colmap, {''}, bw2);
    % xl = [-1 1]*.5;
    for i=1:2
        % set(h(i), 'xlim', xl);
        title(h(i), 'Model');
        ylabel(h(i), ylabel_str{i});
    end
    
    legend(h(1), colstrs, 'location', 'north', 'box', 'off', 'fontsize', fsy);
    
end

if fig_no == 2
    if nargin< 2
        close all;    
        nr = 1;
        nc = 1;
        fsiz = [0 0 .25 .3];
        subplots = 1;
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end    
    h(1) = plot_bar(nr,nc,subplots(1), {mx(4,:)}, {ex(4,:)}, colstrs, {y_str}, [], colmap, {''}, bw2);
    title(h(1), title_str{3});        
end


for i=1:length(h)
    ar = get(h(i),'PlotBoxAspectRatio');
    ar(2) = .9*ar(2);
    set(h(i), 'PlotBoxAspectRatio', ar);
end

end

function [b_vol, b_sto, b_difr, x_labels, y_labels] = regression(outcome, val, vol, sto)

delta = outcome - val;
delta = delta(1:end-1,:);
dm = val(2:end,:) - val(1:end-1,:); % action(t) is an index of m(t-1) after the update

dvol = vol(2:end, :) - vol(1:end-1,:);
dsto = sto(2:end, :) - sto(1:end-1,:);
difr = log(vol(2:end, :)) - log(sto(2:end, :));

X = [];
Z = [];
for i=1:size(delta, 2)
    X = [X; [abs(delta(:, i)) abs(dm(:, i))]];
%     X = [X; [abs(delta(:, i)) abs(dm(:, i)), delta(:,i), dm(:,i)]];
    Z = [Z; [dvol(:, i), dsto(:, i), difr(:, i)]];
end

X = X*[1 1;1 -1]';

x_labels = {'intercept','|PE|+|PU|', '|PE|-|PU|'};
y_labels = {'vol', 'sto', 'difr'};

x = X;
for i=1:size(Z, 2)
    y = Z(:, i);
    b(i, :) = glmfit(x, y);
end

b_vol = b(1, :);
b_sto = b(2, :);
b_difr = b(3, :);

end
