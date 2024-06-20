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
        
            [b_assignment{i}(n, :)] = tools_assignment(dynamics{n}.val, dat.bag);    
    
            [bb] = tools_clustering(dynamics{n}.val, dat.bag);    
            a_clustering{i}(n, :) = bb;        

            [bb] = tools_modulation(dynamics{n}.val, dat.bag);    
            b_modulation{i}(n, :) = bb(5:6);            
        end
    end
    save(fsig, 'lr', 'a_clustering', 'b_assignment','b_modulation', 'models');
end
f = load(fsig);
models = f.models;
b_assignment = f.b_assignment;
b_modulation = f.b_modulation;
a_clustering = f.a_clustering;
lr = f.lr;
for i=1:length(models)
    b_models(:, i) = b_assignment{i}(:, 1);
    mb_models(i, :) = mean(b_assignment{i}(:, 1));
    sb_models(i, :) = serr(b_assignment{i}(:, 1));

    mm_models(i, :) = mean(b_modulation{i});
    sm_models(i, :) = serr(b_modulation{i});
    
    ma_models(i, :) = mean(a_clustering{i});
    sa_models(i, :) = serr(a_clustering{i});    
    mlr_models(i, :) = mean(lr{i});
    slr_models(i, :) = serr(lr{i});

%     [~, pdb(i, :),~, tstat(i)] = ttest(b_models(:, i));
end

fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', 'hpl_assignment'));
f = load(fsig);
b_data = f.b_data(:, 1);
mb_data = mean(f.b_data(:, 1));
sb_data = serr(f.b_data(:, 1));
mx_data = mean(f.x_data);
sx_data = serr(f.x_data);

fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', 'hpl_clustering'));
f = load(fsig);
ma_data = mean(f.a_data);
sa_data = serr(f.a_data);


fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', 'hpl_modulation'));
f = load(fsig);
mm_data = mean(f.b_data(:,5:6));
sm_data = serr(f.b_data(:,5:6));

[~, lr] = model_neutral(experiment);
mlr_data = mean(lr);
slr_data = serr(lr);


mlr = [mlr_data; mlr_models];
slr = [slr_data; slr_models];

ma = [ma_data; ma_models];
sa = [sa_data; sa_models];

mb = [mb_data; mb_models];
sb = [sb_data; sb_models];

mm = [mm_data; mm_models];
sm = [sm_data; sm_models];


% mx = [mx_data; mx_models];
% sx = [sx_data; sx_models];

%--------------------------------------------------------------------------
for i=1:size(mlr,1)
    mx1{i} = [mlr(i,[1 3])' mlr(i,[2 4])']';
    sx1{i} = [slr(i,[1 3])' slr(i,[2 4])']';

    mx2{i} = ma(i, :);
    sx2{i} = sa(i, :);

%     mx3{i} = mx(i, :)';
%     sx3{i} = sx(i, :)';   
    mx3{i} = mb(i, :)';
    sx3{i} = sb(i, :)';       

    mx4{i} = mm(i, :);
    sx4{i} = sm(i, :);
end
x_dots = [b_data, b_models];



%--------------------------------------------------------------------------
if nargin< 2
    close all;
%     nr = 2;
%     nc = 2;
%     fsiz = [0 0 .55 .65];
%     subplots = 1:4;    

    nr = 2;
    nc = 4;
    fsiz = [0 0 .8 .5];  
    subplots = {1:4, 5:8, 1:4, 5:8};
    
    fig(1) = figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

K = nc;

% plots 1
fs = def('fs');
fsy = def('fsy');
lgtitle = def('sto');

glbl = {'Small','Large'};
labels = {'Data', 'Delta-bar-delta', 'Weber-like noisy learning', 'Hierarchical Gaussian filter'};
yl = {[.5 .8], [.3 1], [.5 1], [.2 .35]};

h = plot_bar(nr,nc,subplots{1},mx1, sx1, {'Smal','Large'},[{'Learning rate'}, cell(1,3)]);
lg = legend(h(1),glbl,'fontsize',fs,'location','north','box','off');
title(lg,lgtitle,'fontweight','normal');    
for i=1:K
    set(h(i),'ylim', yl{i});
    xlabel(h(i),def('vol'),'fontsize',fs); 
    title(h(i), labels{i}, 'fontsize', fsy);
end
hh = h;
%------------------------------------
colmap = def('col_yg');
bw1 = .05;
% h = plot_bar(nr,nc,subplots{2}, mx2, sx2, {'Small','Large'},repmat({'Changes in learning rate'},1,K), [], colmap, cell(1,4), bw1); 
h = plot_bar(nr,nc,subplots{2}, mx2, sx2, {'Negative','Positive'},[{'Changes in learning rate'}, cell(1,3)], [], colmap, cell(1,4), bw1); 
for i=1:K
    xlabel(h(i),'Sample Autocorrelation', 'fontsize', fs);
end
hh = [hh h];
%------------------------------------
fsiz = [0 0 .8 .4];
fig(2) = figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);   


colmap = def('col_unique');
bw1 = .2;
xstr = {'20%', '40%', '60%', '80%', '100%'};
yls = sprintf('Relationship between\n |LR| changes and |AC|');

% h = plot_bar(nr,nc,subplots{3}, mx3, sx3, xstr,[{'Changes in |LR|'}, cell(1,3)], [], colmap, cell(1,4), bw1); 
% hh = [hh h];
% xlabel('|AC| (binned)');
% ylim([.2 .55]);
% title('Data');

% h = plot_bar(nr,nc,subplots{3}, mx3, sx3, [],[{'Changes in |LR|'}, cell(1,3)], [], colmap, cell(1,4), bw1); 

for i= 1:4
    h(i) = subplot(nr, nc, subplots{3}(i));    
    plot_raincloud(mx3{i}, sx3{i}, x_dots(:, i), experiment);

    if i==1
        ylabel(yls, 'Interpreter','tex');
    end
end
hh = [hh h];

%------------------------------------

colmap = def('col_yg');
bw1 = .05;
h = plot_bar(nr,nc,subplots{4}, mx4, sx4, {'|PE|+|PU|','|PE|-|PU|'},[{'Learning rate modulation'}, cell(1,3)], [], colmap, cell(1,4), bw1); 
set(h, 'fontsize', fs)
hh = [hh h];

set(hh, 'fontsize', fs)

end