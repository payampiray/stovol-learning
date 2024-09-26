function [st, lr] = model_neutral(experiment, nr, nc, subplots, fig_no)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fname,'file')
    [lr, b] = glm_analysis(experiment);
    save(fname, 'lr', 'b');
end
f = load(fname);
lr = f.lr;
b = f.b;
if nargout>1
    st = [];
    return;
end

%--------------------------------------------------------------------------
eff_lr = lr*[1 1 1 1; -1 -1 1 1; -1 1 -1 1; 1 -1 -1 1]'/2;
eff_simple = b(:, 5:end)*[-1 -1 1 1; -1 1 -1 1; 1 -1 -1 1; 1 1 1 1]'/2;
eff = [eff_lr, eff_simple];
mean_eff = mean(eff);
serr_eff = serr(eff);
% ci_eff = confidence_interval(eff);
[~, p_eff, ci_eff, st] = ttest(eff);
% for i=1:size(eff, 2)
%     bf10(i) = bf.ttest(eff(:, i));
% end

tbl_data = [mean_eff; serr_eff; st.tstat; p_eff];
st.table.rows = {'Mean Effect'; 'S.E.M.'; 't-statistics'; 'P-value'};
st.table.columns = {'PE','PE x Sto', 'PE x Vol', 'PE x Sto x Vol', 'Sto', 'Vol', 'Sto x Vol', 'Intercept'};
st.table.data = tbl_data;


st.p = p_eff;
st.ci = ci_eff;
st.labels = {'PE','PE x Sto', 'PE x Vol', 'PE x Sto x Vol', 'Sto', 'Vol', 'Sto x Vol', 'Intercept'};
st.mean = mean_eff;
st.percent_neg = round(100*mean(eff<0));
% st.bf_null = 1./bf10;

% [R,P,RL,RU] = corrcoef(eff_lr(:, 2:3));
[cr, cp] = corr(eff_lr(:, 2:3));
st.corr_stovol_r = cr(:)';
st.corr_stovol_p = cp(:)';

ma = mean(lr);
ea = serr(lr);
ma = [ma([1 3])' ma([2 4])'];
ea = [ea([1 3])' ea([2 4])'];

f = lr*[-1 -1 1 1;-1 1 -1 1]';
mf = mean(f);
ef = serr(f);

bs = f(:,1);
bv = f(:,2);

%--------------------------------------------------------------------------

colstrs = {'Small','Large'};
glbl = {'Small','Large'};
lgtitle = def('sto');
olbl = {'Small true volatility','Large true volatility'};
ii = [1 3;2 4];
eff_labels = {def('sto'), def('vol')};
eff_labels = {sprintf('True\nstochasticity'), sprintf('True\nvolatility')};


col = def('col');
fn = def('fn');
fsy = def('fsy');
fsl = fsy;
yl = [0 90];
col_br = def('col_yg');

if nargin<2
    close all;
    nr = 1;
    nc = 1;
    fsiz = [0 0 .18 .3];
    subplots = 1;    
    fig_no = 1;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

if fig_no == 1
    h = plot_bar(nr,nc,subplots(1),{ma'},{ea'},colstrs,{'Learning rate coefficient'});
    lg = legend(h(1),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle,'fontweight','normal');
    xlabel(h,def('vol'),'fontsize',fsy); 
    set(h,'ylim',[.5 0.8]);
    set(h, 'ytick', 0.5:.1:.8)
end

% -------------------------------------------------------------------------
% Figure 2
if nargin<1
    nr = 1;
    nc = 1;
    fsiz = [0 0 .40 .25];
    subplots = 1;
    fig_no = 2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

if fig_no == 2

    alf = def('alf');
    fsy = def('fsy');
    
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
    
    for i=1:2
        plot(mf(i)+[-ef(i) ef(i)], yy(i)*ones(1,2),'color','k','linewidth',2); hold on;
    end
    set(gca,'box','off');
    set(gca,'ytick',[]);
    set(gca,'xlim',[-1 1]);
    xlabel('Effect size', 'FontSize', fsy)
    % set(gca, 'yticklabel', eff_labels, 'YTickLabelRotation',-20, 'fontsize',fsy,'TickLabelInterpreter', 'tex')
    
    yy = [0.85, 2.15];
    
    xsA = -1.17;
    for i=1:2
        text(xsA,yy(i),eff_labels{i},'HorizontalAlignment','center','VerticalAlignment','middle','fontsize',fsy,'fontname',fn,'parent',gca);
    end
end

end


function [lr_coeff, b] = glm_analysis(experiment)
[data] = get_data(experiment);

b = nan(length(data), 8);
for n=1:length(data)
    delta_all = []; update_all = [];block_all = [];
    dat = data{n};
    for j=1:4
        
        bucket = dat.bucket(:,j);
        bag = dat.bag(:,j);
        
        update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
        delta = bag - bucket;
        delta = delta(1:end-1);

        update_all = [update_all; update]; %#ok<AGROW> 
        delta_all = blkdiag(delta_all, delta);        
        block_all = blkdiag(block_all,ones(size(update)));
    end
    b(n,:) = glmfit([delta_all block_all],update_all,'normal', 'constant', 'off');
end
lr_coeff = b(:, 1:4);
end