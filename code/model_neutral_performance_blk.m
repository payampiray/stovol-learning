function stats = model_neutral_performance_blk(experiment)
if nargin<1, experiment = 1; end
[data] = get_data(experiment);

e_per_blk = nan(length(data), 4);
for n=1:length(data)
    y = data{n}.bucket(1:end, :);
    x = data{n}.bag(1:end, :);
    e_per_blk(n, :) = mean(abs(x - y).^2);
end
vr_per_blk = var(x);

[q_blk, ~, st_blk] = signrank2(e_per_blk - vr_per_blk);
med_blk = median(e_per_blk - vr_per_blk);

stats.labels = {'block1','block2', 'block3', 'block4'};
stats.median = med_blk;
stats.p = q_blk;
stats.zval = st_blk.zval;
end
