function [stats, figs] = model_neutral_response_time(experiment, nr, nc, subplots, fig_no)
if nargin<1, experiment = 1; end
if nargin<5, fig_no = 1; end
[data] = get_data(experiment);

b_data = nan(length(data), 3);
% x_data = nan(length(data), 10);
for n=1:length(data)
    dat = data{n};                
    [b_data(n, :), x_data(n, :)] = regression_response_time(dat.bucket, dat.bag, dat.response_time);
end

mb = mean(b_data);
sb = serr(b_data);
[~, pb, ci, tstat] = ttest(b_data);
tb = tstat.tstat;

mx = mean(x_data);
sx = serr(x_data);

st = tstat;
st.p = pb;
st.ci = ci;
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


% -------------------------------------------------------------------------
fsy= def('fsy');
fs = def('fs');
yls = sprintf('Relationship between\n response time and |AC|');
xl = 'Relationship between response time and |AC|';
yl = [-1.1 1.1]*10^-4;
colmap = def('col_unique');

if fig_no == 2
    if nargin< 2
        close all;
        nr = 1;
        nc = 1;
        fsiz = [0 0 .2 .3];  
        subplots = 1;    
        
        figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
    end
    figs(1) = gcf;

    
    yl_rt = [475 485]/1000;
    yticklabel = {'0.475', '0.480', '0.485'};
    
    
    subplot(nr, nc, subplots(1));
    plot_raincloud(mb(1), eb(1), b_dots(:, 1), experiment);
    title('Response Time', 'fontsize', fsy);
    ylabel(yls, 'fontsize', fsy);
    ylim(yl);
    
    return;
end
% -------------------------------------------------------------------------
if nargin< 2
    close all;
    nr = 1;
    nc = 1;
    subplots = [1, 2];
end

bw1 = 0.08;

i = 1;
title_str = {'Response time'};

fsiz = [0 0 .16 .22];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);     
plot_bar(nr, nc, subplots(i), {mb(i)}, {eb(i)}, {''}, {yls}, [], colmap, {''}, bw1);
% title('Response Time', 'fontsize', fsy);        
title(title_str{i}, 'fontsize', fsy); 
figs(1) = gcf;

fsiz = [0 0 .3 .22];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);  
subplot(nr, nc, subplots(i));    

config.add_dots = 1;
config.is_right = 1;
config.dens_ratio = 2.5;
config.face_color = colmap;
config.patch_alpha = .4;
config.func_summary = 'median';
y_points = [0.5 3.5];
raincloud1xN_horizontal(b_dots(:, i), y_points, config)    
set(gca, 'fontsize', fs, 'ytick', []);

xlabel(xl, 'fontsize', fsy);
title(title_str{i}, 'fontsize', fsy);    
figs(2) = gcf;


fsiz = [0 0 .2 .22];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);  
wbin = 100/length(mx);
xperc = wbin:wbin:(100);
y_str = 'Response time (s)';
yl_rt = [476 484]/1000;

subplot(nr, nc, subplots(i));
shadedErrorBar(xperc, mx(i, :), sx(i, :), 'lineProps', {'','color', colmap, 'markerfacecolor', [0 0 0]}); hold on;
plot(xperc, mx(i, :), '.','markerfacecolor', [0 0 0], 'markersize', 10); hold on;


set(gca, 'fontsize', fs);
ylabel(y_str, 'fontsize', fsy);
xlabel('|AC|% (binned)', 'fontsize', fsy);
ylim(yl_rt);
title(title_str{i}, 'fontsize', fsy);    
figs(3) = gcf;


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
xt = prctile(acov_abs, 10:10:90);
% xt = prctile(acov_abs, 5:5:95);
xt = [0 xt max(acov_abs)];
y = rt_all;
for i=2:(length(xt))
    t = (acov_abs>=xt(i-1)) & (acov_abs<xt(i));
    x(i-1) = [nanmean(y(t))];
end

end
