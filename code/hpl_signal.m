function stats = hpl_signal(experiment, nr, nc, subplots, fig_no)

if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end

f = get_signals(experiment);
vol = f.vol;
sto = f.sto;
lr = f.lr;

func1 = @mean;
func2 = @serr;
mv = func1(vol, 3);
ms = func1(sto, 3);
ev = func2(vol, 3);
es = func2(sto, 3);

func = @median;
ma1 = func(lr,1);
mv1 = func(vol,1);
ms1 = func(sto,1);
siz = size(ma1);
ma1 = reshape(ma1, siz(2:3))';
mv1 = reshape(mv1, siz(2:3))';
ms1 = reshape(ms1, siz(2:3))';

ma1f = ma1*[-1 -1 1 1;-1 1 -1 1]';
mv1f = mv1*[-1 -1 1 1;-1 1 -1 1]';
ms1f = ms1*[-1 -1 1 1;-1 1 -1 1]';
% [~, qa, ci, stats_lr] = ttest(ma1f);
% stats_lr.ci = ci;
[qa, ~, stats_lr] = signrank2(mv1f);
[qs, ~, stats_sto] = signrank2(ms1f);
[qv, ~, stats_vol] = signrank2(mv1f);
% [qv_diff, ~, wilcox_vol_diff] = signrank(mv1f(:,1)-mv1f(:,2));

labels = {'lr', 'sto', 'vol'};
qq = {qa, qs, qv};
stats_all = {stats_lr, stats_sto, stats_vol};
zval = {stats_sto.zval, stats_vol.zval};
med = {median(ma1f), median(ms1f), median(mv1f)};
m = {mean(ma1f), mean(ms1f), mean(mv1f)};
for i=1:3
    stats.(labels{i}) = stats_all{i};    
    stats.(labels{i}).mean = m{i};
    stats.(labels{i}).med = med{i};
    stats.(labels{i}).p = qq{i};
    stats.(labels{i}).labels = {'Sto Effect', 'Vol Effect'};
end

% [~, tru_sto, tru_vol] = get_data(data_set);
% mds = mean(ms1 - tru_sto);
% mdv = mean(mv1 - tru_vol);
% mdr = mean(mv1./ms1 - tru_vol./tru_sto);
% [~, pds] = ttest(ms1 - tru_sto);
% [~, pdv] = ttest(mv1 - tru_vol);
% [~, pdr] = ttest(mv1./ms1 - tru_vol./tru_sto);
% mds = median(ms1 - tru_sto);
% mdv = median(mv1 - tru_vol);
% [qds] = signrank2(ms1 - tru_sto);
% [qdv] = signrank2(mv1 - tru_vol);

func1 = @mean;
func2 = @(x)serr(x);
mma = func1(ma1);
ema = func2(ma1);
mmv = func1(mv1);
emv = func2(mv1);
mms = func1(ms1);
ems = func2(ms1);


f = ma1*[-1 -1 1 1;-1 1 -1 1]';
mf = median(f);
ef = serr(f);

bs = f(:,1);
bv = f(:,2);

mma = [mma([1 3])' mma([2 4])'];
ema = [ema([1 3])' ema([2 4])'];
mmv = [mmv([1 3])' mmv([2 4])'];
emv = [emv([1 3])' emv([2 4])'];
mms = [mms([1 3])' mms([2 4])'];
ems = [ems([1 3])' ems([2 4])'];

%--------------------------------------------------------------------------
if fig_no == 1
    if nargin< 2
        close all;    
        nr = 1;
        nc = 4;
        fsiz = [0 0 .75 .25];
        subplots = 2:4;
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end   
    
    fsy = def('fsy');
    lgtitle = def('sto');
    
    glbl = {'Small','Large'};
    labels = {'Learning rate', 'Stochasticity estimate', 'Volatility estimate'};
    
% %     h = plot_bar(nr,nc,subplots(1:3),{mma', mms', mmv'}, {ema', ems', emv'}, {'Smal','Large'},labels);
% %     lg = legend(h(1),glbl,'fontsize',fsy,'location','northwest','box','off');
% %     title(lg,lgtitle,'fontweight','normal');    
% %     set(h(1),'ylim',[.4 .8]);
% %     set(h(1), 'ytick', .4:.1:.8);
% %     for i=1:3
% %         xlabel(h(i),def('vol'),'fontsize',fsy); 
% %     end

    h = plot_bar(nr,nc,subplots(1),{mma'}, {ema'}, {'Smal','Large'},labels);
    lg = legend(h(1),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle,'fontweight','normal');    
    set(h(1),'ylim',[.4 .8]);
    set(h(1), 'ytick', .4:.1:.8);
    xlabel(h(1),def('vol'),'fontsize',fsy); 

    return;
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if fig_no == 1.5
    if nargin< 2
        close all;    
        nr = 1;
        nc = 4;
        fsiz = [0 0 .75 .25];
        subplots = 3:4;
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end 

    eff_labels = {sprintf('True\nstochasticity'), sprintf('True\nvolatility')};
    col_br = def('col_yg');    

    alf = def('alf');
    fsy = def('fsy');
    fs = def('fs');
    
    h = subplot(nr, nc, subplots);
    
    rng(0);
    % plot(ones(size(bs)), bs,'.'); hold on;
    % plot(2*ones(size(bv)), bv,'.'); hold on;
    for n=1:size(bs,1)
        plot([bs(n), bv(n)], [1 2]+randn*.05, '-o','color',[1 1 1]*.8,'markersize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]); hold on;
    end
    % yl = [.6 2.4];
    % ylim(yl);
    % plot([0 0],yl,'color','k'); hold on;
    
    wb = .2;
    
    yy = [0.7 2.3];
    for i=1:2
        hb = barh(yy(i), mf(i), wb,'FaceColor',col_br(i,:),'EdgeColor','k','linewidth',1);
        hb.FaceColor = col_br(i, :);
        alpha(alf);
        hold on;
    end
    
%     for i=1:2
%         plot(mf(i)+[-ef(i) ef(i)], yy(i)*ones(1,2),'color','k','linewidth',2); hold on;
%     end
%     set(gca,'xlim',[-1 1.5]);

    set(gca,'box','off');
    set(gca,'ytick',[]);
    xlabel(h, 'Effect size', 'FontSize', fsy)
    % set(gca, 'yticklabel', eff_labels, 'YTickLabelRotation',-20, 'fontsize',fsy,'TickLabelInterpreter', 'tex')
    
    yy = [0.85, 2.15];
    
    xsA = -.77;
    for i=1:2
        text(xsA,yy(i),eff_labels{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',fs,'parent',h);
    end

    return;
end

if fig_no == 0.5
    if nargin< 2
        close all;    
        nr = 1;
        nc = 2;
        fsiz = [0 0 .4 .25];
        subplots = 1:2;
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end   
    
    fsy = def('fsy');
    lgtitle = def('sto');
    
    glbl = {'Small','Large'};
    labels = {'Stochasticity estimate', 'Volatility estimate'};
    
    h = plot_bar(nr,nc,subplots(1:2),{mms', mmv'}, {ems', emv'}, {'Smal','Large'},labels);
    lg = legend(h(2),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle,'fontweight','normal');    
%     set(h(1),'ylim',[.4 .8]);
%     set(h(1), 'ytick', .4:.1:.8);
    for i=1:length(h)
        xlabel(h(i),def('vol'),'fontsize',fsy); 
    end

    return;
end

%--------------------------------------------------------------------------
if nargin< 2
    close all;    
    nr = 2;
    nc = 2;
    fsiz = [0 0 .75 .55];
    subplots = 1:4;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end    
fsy = def('fsy');
fs = def('fs');
cols = def('col');

ii = {[1 3], [2 4]};
label = 'Stochasticity estimate';
lg_labels = {'Small', 'Large'};
lg_title = def('sto');
fig_titles = {sprintf('Small %s', lower(def('vol'))), sprintf('Large %s', lower(def('vol')))};
yl = [0 120];
for i=1:2
    subplot(nr, nc, subplots(i));
    for j=1:2
        hm(j) = plot_ci(ms(:, ii{i}(j)), es(:, ii{i}(j)), cols(j, :)); hold on;
        ylim(yl);
    end
    ylabel(label, 'fontsize', fsy);
    xlabel('Trial', 'fontsize', fsy);
    title(fig_titles{i}, 'fontsize', fsy, 'fontweight', 'normal');
    set(gca, 'box', 'off');

    if i==1
        lg = legend(hm, lg_labels, 'box', 'off', 'location', 'northwest', 'fontsize', fs);
        title(lg, lg_title,'fontweight','normal', 'fontsize', fsy);
    end
end

ii = {[1 2], [3 4]};  
cols = def('col_bp');
label = 'Volatility estimate';
lg_labels = {'Small', 'Large'};
lg_title = def('vol');
fig_titles = {sprintf('Small %s', lower(def('sto'))), sprintf('Large %s', lower(def('sto')))};
yl = [0 180];
for i=1:2
    subplot(nr, nc, subplots(2+i));
    for j=1:2
        hm(j) = plot_ci(mv(:, ii{i}(j)), ev(:, ii{i}(j)), cols(j, :)); hold on;
        ylim(yl);
    end
    ylabel(label, 'fontsize', fsy);
    xlabel('Trial', 'fontsize', fsy);
    title(fig_titles{i}, 'fontsize', fsy, 'fontweight', 'normal');
    set(gca, 'box', 'off');
    
    if i==1
        lg = legend(hm, lg_labels, 'box', 'off', 'location', 'northwest', 'fontsize', fs);
        title(lg, lg_title,'fontweight','normal', 'fontsize', fsy);
    end
end



end

function signals = get_signals(data_set)

fname = fullfile('..',sprintf('experiment%d', data_set), 'model_hpl.mat');
f = load(fname);
N = length(f.dynamics);
for n=1:N
    lr(:,:,n) = f.dynamics{n}.lr;
    vol(:,:,n) = f.dynamics{n}.vol;
    sto(:,:,n) = f.dynamics{n}.sto;
end
signals = struct('vol', vol, 'sto', sto, 'lr', lr);
end
