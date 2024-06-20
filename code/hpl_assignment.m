function stats = hpl_assignment(experiment, nr, nc, subplots, fig_no)

if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end
[data] = get_data(experiment);

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
fsig = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fsig, 'file')
    for n=1:length(data)
        dat = data{n};                
        [b_data(n, :), x_data(n,:)] = tools_assignment(dat.bucket, dat.bag);
    
        dynamics = f.dynamics{n};
        [b_model(n, :), x_model(n, :), bv2s(n, :)] = tools_assignment(dynamics.val, dat.bag, dynamics.vol, dynamics.sto);
    
        save(fsig, 'b_data', 'x_data', 'b_model', 'x_model', 'bv2s');
    end
end

f = load(fsig);
b_data = f.b_data;
x_data = f.x_data;
b_model = f.b_model;
x_model = f.x_model;
bv2s = f.bv2s;

b = {b_model, b_data, bv2s};
for i=1:length(b)
    xdb(i, :) = median(b{i});
    mdb(i, :) = mean(b{i});
    sdb(i, :) = serr(b{i});
    [qdb(i, :),~, wilcox(i)] = signrank2(b{i});
    [~, pdb(i, :),~, tstat(i)] = ttest(b{i});
    tval(i, :) = tstat(i).tstat;
    se(i, :) = tstat(i).sd;
end

x = {x_model, x_data};
for i=1:length(x)
    mx(i, :) = mean(x{i});
    sx(i, :) = serr(x{i});
end



y_labels = {'model', 'data', 'b2v'};
for i=1:3
    st = tstat(i);
    st.p = pdb(i, :);
    st.mean = mdb(i, :);
    st.sem = sdb(i, :);
    st.columns = {'|AC|', 'AC', 'intercept'};
    stats.(y_labels{i}) = st;
end

% tbl_data = [mdb(2,:); sdb(2,:); tval(2,:); pdb(2,:)];
% 
% stats.table(1).data = tbl_data;
% stats.table(1).columns = {'|AC|', 'AC', 'intercept'};
% stats.table(1).rows = {'Estimate'; 'SE'; 't-value'; 'P-value'};
% 
% stats.table(2).data = [mx; sx];
% stats.table(2).columns = {'20%', '40%', '60%', '80%', '100%'};
% stats.table(2).rows = {'Mean'; 'SE'};
%--------------------------------------------------------------------------
mb = mdb(1:3, 1);
eb = sdb(1:3, 1);

mx = mx(2, :);
ex = sx(2, :);
b_dots = [b_model(:, 1), b_data(:, 1), bv2s(:, 1)];

if nargin< 2
    close all;
    nr = 3;
    nc = 2;
    fsiz = [0 0 .4 .9];  
    subplots = [1 3 2 4];
    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

fsy = def('fs');
colmap = def('col_unique');
bw1 = 0.04;
yl = [-19 9]*10^-4;

yls{1} = sprintf('Relationship between\n |LR| changes and |AC|');
yls{2} = sprintf('Relationship between\n |LR| changes and |AC|');
yls{3} = sprintf('Relationship between\nlog %s and |AC|', '\DeltaV/\DeltaS');
title_str = {'Model', 'Data', 'Model'};

ii = 1:3;
if fig_no == 0.5, ii = 1:2; end


for i=ii
    subplot(nr, nc, subplots(i));    
    plot_raincloud(mb(i), eb(i), b_dots(:, i), experiment);
    ylabel(yls{i}, 'Interpreter','tex');
    title(title_str{i});
    
    if i~=3
        ylim(yl);
    end
end

if fig_no == 0.5
    return;
end

cols = repmat(colmap(1, :), 10, 1);
% xstr = {'10%', '20', '30', '40', '50', '60', '70', '80', '90', '100'};
xstr = {'10%', '', '30%', '', '50%', '', '70%', '', '90%', ''};
xstr = {'20%', '40%', '60%', '80%', '100%'};
y_str = 'Changes in |LR|';
plot_bar(nr, nc, subplots(4), {mx}, {ex}, '', {y_str}, [], cols, {''}, bw1, 1);
set(gca, 'XTickLabel', xstr, 'fontsize', fsy);
xlabel('|AC| (binned)');
ylim([.2 .55]);
title('Data');

end
