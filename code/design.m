function h = design(nr, nc, subplots, fig_no)

[data, true_sto, true_vol, ~, bird] = get_data(1);

state = bird;
outcome = data{1}.bag;


[lr] = kalman(true_vol, true_sto, outcome);


ma = mean(lr);
ea = ma*0;
% ea = serr(mlr);
ma = [ma([1 3])' ma([2 4])'];
ea = [ea([1 3])' ea([2 4])'];


v = var(outcome);
x = outcome(1:end-1, :);
y = outcome(2:end, :);
a = corr(x,y);
a = diag(a)';
v = [v([1 3])' v([2 4])'];
a = [a([1 3])' a([2 4])'];

%--------------------------------------------------------------------------
if nargin< 2
    close all;    
    nr = 2;
    nc = 1;
    fsiz = [0 0 .4 .4];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 1;
end

colstrs = {'Small','Large'};
glbl = {'Small','Large'};
lgtitle = def('sto');
olbl = {'Small true volatility','Large true volatility'};
ii = [1 3;2 4];


col = def('col');
fn = def('fn');
fsy = def('fsy');
fsl = fsy;
yl = [0 90];

h = nan(1,length(subplots));
if fig_no == 1
    for i=1:2
        h(i) = subplot(nr,nc,subplots(i));
        hp = nan(1,2);            
        plot(state(:,i),'-','linewidth',1,'color','k');

        hold on;
        for j=1:2
            ij = ii(i,j);
            hp(j) = plot(outcome(:,ij),'-','color',col(j,:),'linewidth',1); hold on;
            plot(outcome(:,ij),'.','color',col(j,:),'markersize',10); hold on;
        end    
        set(gca,'fontname',fn,'box','off');
        title(olbl{i},'fontsize',fsy,'fontweight','normal');     
        if i==1
            lg = legend(hp,glbl,'fontsize',fsy,'box','off','orientation','horizontal','location','southeast');    
            lgt = title(lg,lgtitle,'fontweight','normal');
        end
        ylim(yl);        
        ylabel('Observations','fontsize',fsy);
        if i==2
            xlabel('Trial','fontsize',fsy);
        end
    end
end

%--------------------------------------------------------------------------
if nargin<1
    nr = 1;
    nc = 2;
    fsiz = [0 0 .4 .25];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 2;
end

if fig_no == 2
    h = plot_bar(nr,nc,subplots(1:2),{v',a'},{0*v',0*a'},colstrs,{'Variance','Autocorrelation'});
    lg = legend(h(1),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle,'fontweight','normal');
    for i=1:2
        xlabel(h(i),def('vol'),'fontsize',fsy); 
    end
    set(h(2),'ylim',[0 1]);
end

%--------------------------------------------------------------------------
if nargin<1
    nr = 1;
    nc = 1;
    fsiz = [0 0 .18 .25];
    subplots = 1;        
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
    fig_no = 3;
end

if fig_no == 3
    colstrs = {'Small','Large'};
    glbl = {'Small','Large'};
    lgtitle = def('sto');
    fsy = def('fsy');
    
    h = plot_bar(nr,nc,subplots(1),{ma'},{ea'},colstrs,{'Optimal learning rate'});
    lg = legend(h(1),glbl,'fontsize',fsy,'location','northwest','box','off');
    title(lg,lgtitle,'fontweight','normal');
    for i=1
        xlabel(h(i),def('vol'),'fontsize',fsy); 
    end
    % set(h(1),'ylim',[0 .8]);

end

end


function [lr] = kalman(true_vol, true_sto, outcome)

lambda = true_vol./true_sto;

ndim = size(outcome,2);
N = size(outcome,1);
lr = nan(N,ndim);

m = 60+zeros(1,ndim);
r = 1*ones(1,ndim);

for t=1:N

    a = (r+lambda)./(r+1+lambda);
    m = m + a.*(outcome(t, :)-m);
    r = (1-a).*(r+lambda); 

    lr(t, :) = a;
end
end

