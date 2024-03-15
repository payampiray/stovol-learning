function stats = hpl_fluctuations(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('model_neutral_fluctuations.mat'));
f = load(fname);
xb = f.x;

fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fname, 'file')
    [data] = get_data_model(experiment);
    [x, e] = tools_fluctuations(data);
    save(fname, 'x', 'e');
end
f = load(fname);
mb = f.x;
% me = f.e;

[c, q] = corr(xb, mb, 'type', 'spearman');
[c0, q0] = corr(xb(:), mb(:), 'type', 'spearman');
c = diag(c)';
q = diag(q)';

stats.all.corr_spearman_c = c0;
stats.all.corr_spearman_p = q0;
stats.blk.corr_spearman_c = c;
stats.blk.corr_spearman_p = q;

tbl.rows = {'blk1','blk2','blk3','blk4'};
tbl.columns = {'correlation','pval'};
tbl.data = [c' q'];

stats.table = tbl;
end

function [data_model] = get_data_model(experiment)
f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
dynamics = f.dynamics;

data = get_data(experiment);
data_model = cell(size(data));
for n=1:length(data_model)
    val = dynamics{n}.val;
    outcome = data{n}.bag;
    data_model{n} = struct('bucket', val, 'bag', outcome);        
end    

end


