function tbl = other_weber_fit(experiment)
if nargin<1, experiment = 1; end


[lme, parameters] = get_weber(experiment);
x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'alpha', 'zeta', 'sigma'};
tbl.columns = {'25%', '50%', '75%'};

bmc = other_bmc_pf(lme);
end

function [lme, parameters] = get_weber(experiment)

fname = fullfile('..',sprintf('experiment%d', experiment), 'model_weber.mat');
f = load(fname);
parameters = f.parameters.Variables;
num_parameters = width(f.parameters);
logf = -.5*f.loglik -0.5*log(2*pi);
num_observations = 200;

lme = logf -0.5*num_parameters*log(num_observations);
end
