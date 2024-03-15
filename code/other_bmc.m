function st = other_bmc(experiment)
if nargin<1, experiment = 1; end

models = {'model_rw', 'model_dbd1', 'model_dbd4','model_weber', 'model_hgf', 'model_kfem'};

fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hpl.mat');
f = load(fname);
lme_pf = f.lme;

for i=1:length(models)
    fname = fullfile('..',sprintf('experiment%d', experiment), sprintf('%s.mat', models{i}));
    f = load(fname);
    tbl_model = bmc_run(lme_pf, f.lme);
    if i>1
        tbl.data = [tbl.data; tbl_model.data];
    else
        tbl = tbl_model;
    end
end
tbl.rows = models;

st.table = tbl;
end

function tbl = bmc_run(lme_pf, other_lme)


lme = [other_lme lme_pf];

[~, mf, ~, pxp, ~] = cbm_spm_BMS(lme);

x = [mf pxp(2)];

tbl.data = x;
tbl.columns = {'MF:other', 'MF:pf', 'PXP:pf'};
end
