function stats = hpl_bmc(experiment)
if nargin<1, experiment = 1; end


fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hpl.mat');
f = load(fname);
[lme(:, 1)] = f.lme;

fname = fullfile('..',sprintf('experiment%d', experiment),'model_kf4.mat');
f = load(fname);
[lme(:, 2)] = f.lme;

[~, mf, ~, pxp, ~] = cbm_spm_BMS(lme);

model_names = {'pf', 'kf'};
stats = struct('model_frequency', mf, 'pxp', pxp, 'labels', {model_names});

end

