function hh = other_signals(experiment)
if nargin<1, experiment = 1; end

data = get_data(experiment);

models = {'model_dbd1', 'model_weber', 'model_hgf'};
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));

if ~exist(fsig, 'file')
    for i=1:length(models)
        fname = fullfile('..',sprintf('experiment%d', experiment), sprintf('%s.mat', models{i}));
        f = load(fname);
        dynamics = f.dynamics;
        for n=1:length(dynamics)
            dat = data{n};
            lr{i}(n, :) = mean(dynamics{n}.lr);
        
            [bb] = tools_modulation(dynamics{n}.val, dat.bag);    
            b_modulation{i}(n, :) = bb(5:6);
    
            [bb] = tools_clustering(dynamics{n}.val, dat.bag);    
            a_clustering{i}(n, :) = bb;        
        end
    end
    save(fsig, 'lr', 'b_modulation','a_clustering', 'models');
end
f = load(fsig);
models = f.models;
b_modulation = f.b_modulation;
a_clustering = f.a_clustering;
lr = f.lr;
for i=1:length(models)
    mb_models(i, :) = mean(b_modulation{i});
    sb_models(i, :) = serr(b_modulation{i});
    ma_models(i, :) = mean(a_clustering{i});
    sa_models(i, :) = serr(a_clustering{i});    
    mlr_models(i, :) = mean(lr{i});
    slr_models(i, :) = serr(lr{i});
end

fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', 'hpl_modulation'));
f = load(fsig);
mb_data = mean(f.b_data(:,5:6));
sb_data = serr(f.b_data(:,5:6));

fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', 'hpl_clustering'));
f = load(fsig);
ma_data = mean(f.a_data);
sa_data = serr(f.a_data);

[~, lr] = model_neutral(experiment);
mlr_data = mean(lr);
slr_data = serr(lr);


mlr = [mlr_data; mlr_models];
slr = [slr_data; slr_models];

ma = [ma_data; ma_models];
sa = [sa_data; sa_models];

mb = [mb_data; mb_models];
sb = [sb_data; sb_models];

%--------------------------------------------------------------------------
for i=1:size(mlr,1)
    mx1{i} = [mlr(i,[1 3])' mlr(i,[2 4])']';
    sx1{i} = [slr(i,[1 3])' slr(i,[2 4])']';

    mx2{i} = ma(i, :);
    sx2{i} = sa(i, :);

    mx3{i} = mb(i, :);
    sx3{i} = sb(i, :);    
end



%--------------------------------------------------------------------------
if nargin< 2
    close all;
%     nr = 2;
%     nc = 2;
%     fsiz = [0 0 .55 .65];
%     subplots = 1:4;    

    nr = 3;
    nc = 4;
    fsiz = [0 0 .8 .9];  
    subplots = {1:4, 5:8, 9:12};
    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

K = nc;

% plots 1
fs = def('fs');
fsy = def('fsy');
lgtitle = def('sto');

glbl = {'Small','Large'};
labels = {'Data', 'Delta-bar-delta', 'Weber-like noisy learning', 'Hierarchical Gaussian filter'};
yl = {[.5 .8], [.6 1], [.5 .64], [.2 .35]};

h = plot_bar(nr,nc,subplots{1},mx1, sx1, {'Smal','Large'},[{'Learning rate'}, cell(1,3)]);
lg = legend(h(1),glbl,'fontsize',fs,'location','north','box','off');
title(lg,lgtitle,'fontweight','normal');    
for i=1:K
    set(h(i),'ylim', yl{i});
    xlabel(h(i),def('vol'),'fontsize',fsy); 
    title(h(i), labels{i}, 'fontsize', fsy);
end
hh = h;
%------------------------------------
colmap = def('col_yg');
bw1 = .05;
% h = plot_bar(nr,nc,subplots{2}, mx2, sx2, {'Small','Large'},repmat({'Changes in learning rate'},1,K), [], colmap, cell(1,4), bw1); 
h = plot_bar(nr,nc,subplots{2}, mx2, sx2, {'Small','Large'},[{'Changes in learning rate'}, cell(1,3)], [], colmap, cell(1,4), bw1); 
for i=1:K
    xlabel(h(i),'Sample Autocorrelation', 'fontsize', fsy);
end
hh = [hh h];
%------------------------------------
colmap = def('col_yg');
bw1 = .05;
h = plot_bar(nr,nc,subplots{3}, mx3, sx3, {'|PE|+|PU|','|PE|-|PU|'},[{'Learning rate modulation'}, cell(1,3)], [], colmap, cell(1,4), bw1); 
hh = [hh h];
%------------------------------------


end