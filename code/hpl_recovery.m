function st = hpl_recovery(nr, nc, subplots)

nn = 1:100;
fname = fullfile('..',sprintf('experiment%d', 1), 'hpl_sim.mat');   
f = load(fname);
parameters = table2array(f.fitted_parameters(nn, :));
true_parameters = table2array(f.true_parameters(nn, :));

N = [sign(parameters(:,1)-parameters(:,2)) sign(true_parameters(:,1)-true_parameters(:,2))];
idx = (N(:,1)~=N(:,2));
m_bad = mean(idx);

r = corr(parameters,true_parameters,'type','spearman');
r = diag(r);

tbl_data = r';
st.table.rows = {'Correlation'};
st.table.columns = {'$\mu_s$', '$\mu_v$', '$\sigma^2$'};
st.table.data = tbl_data;

%==========================================================================
params = true_parameters-parameters;
% mo = mean(true);

if nargin<1
    close all;    
    nr = 1;
    nc = 3;
    fsiz = [0.1    0.0800    .75    .3];
    subplots = 1:3;
%     nc = 4;
%     fsiz = [0.1    0.0800    .95    .3];
%     subplots = 1:4;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end


pnames = {'$\mu_s$', '$\mu_v$', '$\sigma^2$'};
% labels = {'True', 'Fitted'};

fsy = def('fsy');
fsx = def('fsy')+4;
% bandwidth = 0.2;

col_curve = [.5 .2 .1];
col_line = [.8 .2 .2];
alf = .2;

for i=1:nc

    h(i) = subplot(nr,nc,subplots(i));
    [fq,xq, bw] = ksdensity(params(:,i) ,'function','pdf','kernel','normal');
    plot(xq,fq, 'color', col_line, 'linewidth', 1); hold on;
    % [fq,xq] = ksdensity(q_pmmh(burnin+1:end));
    % plot(xq,fq,'r-');
%     ylim([0 1.1]);
%     xlim([-1 2]);

    yl=get(gca,'ylim');
    xlabel('True – Fitted','fontsize',fsy);
    ylabel('Empirical distribution','fontsize',fsy);
    ht = title(pnames{i},'Interpreter','latex','fontsize',fsx);
    
%     mxq = mean(params(:,i));
%     sxq = std(params(:,i))*[-1 1];
%     mxq = median(params(:,i));
%     sxq = se_median(params(:,i))*[-1 1];

    mxq = median(params(:,i));
    sxq = quantile(params(:,i), [.25 .75]) - mxq;    

    [~,t1] = min( abs(xq - (mxq+sxq(1)) ));
    [~,t2] = min( abs(xq - (mxq+sxq(2)) ));

    parameters = xq(t1:t2);
    x2 = [parameters, fliplr(parameters)];
%     yl = get(gca,'ylim');
    inBetween = [yl(1)*ones(1,length(parameters)), yl(2)*ones(1,length(parameters))];
    fill(x2, inBetween, col_curve, 'FaceAlpha', alf, 'EdgeColor', col_curve,'EdgeAlpha', alf); hold on;       
    
    plot([mxq mxq],yl*2, 'color', col_line, 'linewidth',2);
%     plot([0 0],yl*2,'k-','linewidth',1);
    set(gca,'ylim',yl);
    set(gca,'xlim', [-1 1]);
end

end
