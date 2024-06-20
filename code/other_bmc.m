function st = other_bmc(experiment)
if nargin<1, experiment = 1; end

models = {'model_rw', 'model_dbd1', 'model_dbd4','model_weber', 'model_hgf', 'model_kfem'};
models = {'model_rw', 'model_AC'};
fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hpl.mat');
f = load(fname);
lme_pf = f.lme;

for i=1:length(models)
    fname = fullfile('..',sprintf('experiment%d', experiment), sprintf('%s.mat', models{i}));
    f = load(fname);
    lme(:, i) = f.lme;
    tbl_model = bmc_run(lme_pf, f.lme);
    if i>1
        tbl.data = [tbl.data; tbl_model.data];
    else
        tbl = tbl_model;
    end
end

lme = [lme, lme_pf];
[~, mf, ~, pxp, ~] = cbm_spm_BMS(lme);

tbl.rows = [models, 'PF'];

tbl.data = [tbl.data; nan(1, 3)];
tbl.data = [tbl.data mf' pxp'];
tbl.columns = [tbl.columns 'MF', 'PXP'];

st.table = tbl;
end

function tbl = bmc_run(lme_pf, other_lme)


lme = [other_lme lme_pf];

[~, mf, ~, pxp, ~] = cbm_spm_BMS(lme);

x = [mf pxp(2)];

tbl.data = x;
tbl.columns = {'MF:alternative', 'MF:pf', 'PXP:pf'};
end
