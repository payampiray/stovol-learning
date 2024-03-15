function st1 = model_neutral_1st_block(experiment)
if nargin<1, experiment = 1; end

[data, ~, ~, ~, randomized] = get_data(experiment);

ur = unique(randomized(:,1), 'rows');
x = zeros(length(data), size(ur, 1));
for j=1:size(ur, 1)
    n = randomized(:, 1) == ur(j, :);
    x(n, j) = 1;
end

x = x*[-1 -1 1 1;-1 1 -1 1;1 -1 -1 1]';

fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));
if ~exist(fname, 'file')
    glme = run_mixed_glm(data, x);
    save(fname, 'glme');
end
f = load(fname);
st = f.glme;

idx = [2 7 8 10 3 4 9 1 5 6];
x = [st.Coefficients.Estimate st.Coefficients.SE st.Coefficients.tStat st.Coefficients.pValue]';
st1.table.data = x(:, idx);
st1.table.columns = st.CoefficientNames(idx);
st1.table.rows = {'Coeff', 'SEM','tval','pval'};

for i=1:4
    st1.(st1.table.rows{i}) = st1.table.data(i, :);
end
st1.labels = st1.table.columns;

end

function glme = run_mixed_glm(data, session)

cc = [1 1 1 1;-1 -1 1 1;-1 1 -1 1;1 -1 -1 1];

dd =[];
bb =[];
nn = [];
uu = [];
ss = [];
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
    dd = [dd; delta_all*cc'];
    bb = [bb; block_all*cc(2:4,:)'];
    uu = [uu; update_all];
    nn = [nn; n*ones(size(update_all))];
    ss = [ss; session(n,:).*ones(size(update_all))];
end

labels = {'PE','PE x Sto', 'PE x Vol', 'PE x Sto x Vol', 'Sto', 'Vol', 'Sto x Vol', 'update', 'subj', 'S_Sto', 'S_Vol', 'S_Sto_Vol'};
dat = array2table([dd bb uu nn ss], 'VariableNames',labels);

glme = fitglme(dat,'update ~ 1 + PE*Sto*Vol + S_Sto + S_Vol + (PE*Sto*Vol | subj) + (1| subj)');
end

