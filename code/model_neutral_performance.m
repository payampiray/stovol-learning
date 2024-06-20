function stats = model_neutral_performance(experiment)
if nargin<1, experiment = 1; end
[data] = get_data(experiment);

[~, lr] = model_neutral(experiment);
lr_eff = lr*[-1 -1 1 1; -1 1 -1 1]';
fx = lr_eff(:, 1:2).*[-1 1];
fx(:,3) = sum(fx>0, 2)>0;

[err_performance] = nan(size(data));
for n=1:length(data)
    y = data{n}.bucket(1:end, :);
    x = data{n}.bag(1:end, :);
    err_performance(n) = median((mean(abs(x - y)).^2));
end

for i=1:size(fx,2)
    groups = {fx(:,i)>0, (fx(:,i)<=0)};
    [q(i),~, wilxoc(i)] = ranksum(err_performance(groups{2}, :), err_performance(groups{1}, :));
    wilcox_zval(i) = wilxoc(i).zval;
    med_err(i) = median(err_performance(groups{2}, :)) - median(err_performance(groups{1}, :));       
end

adpative = (fx(:,1)>0) & (fx(:, 2)>0);
groups = {adpative, fx(:, 1)<0, fx(:, 2)<0};
for i=1:3
    tbl_data(:, i) = prctile(err_performance(groups{i}), [25 50 75]);   
end

stats.p = q;
stats.zval = wilcox_zval;
stats.labels = {'Stochasticity maladaptive','Volatility maladaptive', 'any maladaptive'};
stats.median = med_err;


stats.table.rows = {'25%'; '50%'; '75%'};
stats.table.columns = {'Adaptive', 'Stochasticity maladaptive','Volatility maladaptive'};
stats.table.data = tbl_data;
end
