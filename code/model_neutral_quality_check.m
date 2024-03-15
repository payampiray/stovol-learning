function model_neutral_quality_check(experiment)
if nargin<1, experiment = 1; end
fname = fullfile('..',sprintf('experiment%d', experiment),sprintf('%s.mat', mfilename));

[data] = get_data(experiment);
if ~exist(fname,'file')
    [outliers, b_PE, p_PE] = quality_check(data);
    save(fname, 'outliers', 'b_PE', 'p_PE');
end
fname = load(fname);
outliers = fname.outliers;
b_PE = fname.b_PE;
p_PE = fname.p_PE;

fprintf('Number of outliers: %d\n', sum(outliers));
fprintf('Number of adverse PE effect: %d\n', sum(b_PE<0));
fprintf('Number of significant positive PE effect: %d\n', sum(p_PE<0.001));

end


function [outliers, b_PE, p_PE] = quality_check(data)
updates = nan(length(data),1);
b_PE = nan(length(data),1);
p_PE = nan(length(data),1);
for i=1:length(data)
    dd = []; dx = [];ff = [];
    dat = data{i};
    for j=1:4
        
        bucket = dat.bucket(:,j);
        bag = dat.bag(:,j);
        
        update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
        delta = bag - bucket;
        delta = delta(1:end-1);

        dx = [dx;update]; %#ok<AGROW> 
        dd = blkdiag(dd,delta);        
        ff = blkdiag(ff,ones(size(update)));        
    end
    cc = [0.25 0.25 0.25 0.25; -1 -1 1 1;-1 1 -1 1;1 -1 -1 1]';
    X1 = dd*cc;
    [bb, ~, st] = glmfit([X1 ff],dx,'normal', 'constant', 'off');
    b_PE(i,:) = bb(1);
    p_PE(i,:) = st.p(1);

    updates(i,:) = mean(abs(dx));
end

mu = mean(updates);
su = std(updates);
outliers = updates<(mu - 3*su);

end
