function stats = model_neutral_response_time(experiment, nr, nc, subplots)
if nargin<1, experiment = 1; end
[data] = get_data(experiment);

b_data = nan(length(data), 3);
x_data = nan(length(data), 5);
for n=1:length(data)
    dat = data{n};                
    [b_data(n, :), x_data(n, :)] = regression_response_time(dat.bucket, dat.bag, dat.response_time);
end

mb = mean(b_data);
sb = serr(b_data);
[~, pb, ~, tstat] = ttest(b_data);
tb = tstat.tstat;

mx = mean(x_data);
sx = serr(x_data);

st = tstat;
st.p = pb;
st.mean = mb;
st.sem = sb;
st.columns = {'|AC|', 'AC', 'intercept'};
stats.data = st;


tbl_data = [mb; sb; tb; pb];
stats.table(1).data = tbl_data;
stats.table(1).columns = {'|AC|', 'AC', 'intercept'};
stats.table(1).rows = {'Estimate'; 'SE'; 't-value'; 'P-value'};

stats.table(2).data = [mx; sx];
stats.table(2).columns = {'20%', '40%', '60%', '80%', '100%'};
stats.table(2).rows = {'Mean'; 'SE'};
% -------------------------------------------------------------------------
mb = mb(1);
eb = sb(1);

b_dots = b_data(:, 1);

if nargin< 2
    close all;
    nr = 1;
    nc = 2;
    fsiz = [0 0 .4 .3];  
    subplots = 1:6;    
    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

colmap = def('col_unique');
bw1 = 0.04;
fs = def('fs');

yl = [-1.1 1.1]*10^-4;
yl_rt = [475 485]/1000;
yticklabel = {'0.475', '0.480', '0.485'};
yls = sprintf('Relationship between\n response time and |AC|');

subplot(nr, nc, subplots(1));
plot_raincloud(mb(1), eb(1), b_dots(:, 1), experiment);
title('Response Time Data');
ylabel(yls);
ylim(yl)

cols = repmat(colmap(1, :), 10, 1);
xstr = {'20%', '40%', '60%', '80%', '100%'};
y_str = 'Response time (s)';
plot_bar(nr, nc, subplots(2), {mx}, {sx}, '', {y_str}, [], cols, {''}, bw1, 1);
set(gca,'XTickLabel', xstr, 'YTickLabel', yticklabel, 'fontsize', fs);
xlabel('|AC| (binned)');
title('Response Time Data');

ylim(yl_rt)
end

function [b, x] = regression_response_time(action, outcome, rt)
acov = [];
rt_all = [];
for j=1:4
    
    bucket = action(:,j);
    bag = outcome(:,j);
    
    delta = bag - bucket;
    delta = delta(1:end-1);

    acov = cat(1, acov, delta.*[0; delta(1:end-1)]);

    rt0 = rt(2:end, j);
    isoutlier_rt = isoutlier(rt0);
    rt0(isoutlier_rt) = nan;
    rt_all = cat(1, rt_all, rt0);
end

rt_all = rt_all/10000; %ms
acov_abs = abs(acov);

b = glmfit([acov_abs acov], rt_all);
b = [b(2:end); b(1)];

xt = prctile(acov_abs, 20:20:80);
xt = [0 xt max(acov_abs)];
y = rt_all;
for i=2:(length(xt))
    t = (acov_abs>=xt(i-1)) & (acov_abs<xt(i));
    x(i-1) = [nanmean(y(t))];
end

end
