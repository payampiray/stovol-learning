function [data, true_sto, true_vol, randomized, bird] = get_data(experiment)
% % % experiment = 1;
fname = fullfile('..',sprintf('experiment%d', experiment), 'data.mat');
f = load(fname);
meta_data = f.meta_data;
data = f.data;
true_sto = meta_data.true_sto;
true_vol = meta_data.true_vol;
bag = meta_data.bag;
bird = meta_data.bird;

data4fit = data;
for n=1:length(data)
    data{n}.bag = bag;
    randomized(n, :) = data{n}.randomization_order;
end    

end