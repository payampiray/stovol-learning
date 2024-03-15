function [b, b_labels] = tools_agnostic_regression(action, outcome)

update_all = [];
delta_all = [];
blocks = [];
PE = [];
PU = [];
for j=1:4
    
    bucket = action(:,j);
    bag = outcome(:,j);
    
    update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
    delta = bag - bucket;
    delta = delta(1:end-1);

    update_all = [update_all; update]; %#ok<AGROW> 
    delta_all = blkdiag(delta_all, delta);        
    blocks = blkdiag(blocks,ones(size(update)));

    PE = blkdiag(PE, [0; abs(delta(1:end-1))]);
    PU = blkdiag(PU, [0; abs(update(1:end-1))]);

end
contrasts = [1 1 1 1; -1 -1 1 1; -1 1 -1 1; 1 -1 -1 1]'/2;

delta_eff = delta_all*contrasts;
PE = PE*[1 1 1 1]';
PU = PU*[1 1 1 1]';

x1 = .5*(PE+PU);
x2 = .5*(PE-PU);

contrasts = [-1 -1 1 1; -1 1 -1 1; 1 -1 -1 1;1 1 1 1]'/2;
block_eff = blocks*contrasts;

[b, ~, st] = glmfit([delta_eff delta_eff(:,1).*x1 delta_eff(:,1).*x2 x1 x2 block_eff],update_all,'normal', 'constant', 'off');
b_labels = {'D', 'D x S', 'D x V', 'D x S x V', 'D x (|PE|+|PU|)', 'D x (|PE|-|PU|)', '|PE|+|PU|', '|PE|-|PU|', 'S', 'V', 'S x V', 'intercept'};
end
