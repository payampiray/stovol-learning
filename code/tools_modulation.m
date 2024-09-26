function [b, b_labels] = tools_modulation(action, outcome)

update_all = [];
delta_all = [];
blocks = [];
PE_abs = [];
PU_abs = [];
lr_pre = [];
for j=1:4
    
    bucket = action(:,j);
    bag = outcome(:,j);
    
    update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
    delta = bag - bucket;
    delta = delta(1:end-1);
    lr = update./delta;

    update_all = [update_all; update]; %#ok<AGROW> 
    delta_all = blkdiag(delta_all, delta);        
    blocks = blkdiag(blocks,ones(size(update)));

    lr_pre = [lr_pre; [0; lr(1:end-1)]];

    PE_abs = blkdiag(PE_abs, [0; abs(delta(1:end-1))]);
    PU_abs = blkdiag(PU_abs, [0; abs(update(1:end-1))]);
end
contrasts = [1 1 1 1; -1 -1 1 1; -1 1 -1 1; 1 -1 -1 1]'/2;

delta_eff = delta_all*contrasts;
PE_abs = PE_abs*[1 1 1 1]';
PU_abs = PU_abs*[1 1 1 1]';

x1 = .5*(PE_abs+PU_abs);
x2 = .5*(PE_abs-PU_abs);

contrasts = [-1 -1 1 1; -1 1 -1 1; 1 -1 -1 1;1 1 1 1]'/2;
block_eff = blocks*contrasts;

[b] = glmfit([delta_eff delta_eff(:,1).*x1 delta_eff(:,1).*x2 delta_eff(:,1).*lr_pre x1 x2 lr_pre block_eff],update_all,'normal', 'constant', 'off');
b_labels = {'D', 'D x S', 'D x V', 'D x S x V', 'D x (|PE|+|PU|)', 'D x (|PE|-|PU|)', 'D x (PU/PE)', '|PE|+|PU|', '|PE|-|PU|', 'PU/PE', 'S', 'V', 'S x V', 'intercept'};
end