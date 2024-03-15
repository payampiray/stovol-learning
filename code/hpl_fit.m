function hpl_fit(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hpl.mat');
    
if ~exist(fname, 'file')
    data = get_data(experiment);
    config = struct('num_bayespoint', 200, 'num_seeds', 5, 'seed', 0, 'num_inits', 10);    
    model = @hpl_model;
    
    [~, vars] = model([], [], 'config');
    opt_var = vars.opt_var;
    InitialX = vars.InitialX;    

    loglik = nan(length(data), 1);
    dynamics = cell(length(data), 1);

    fit_dir = fullfile('..',sprintf('experiment%d', experiment), 'fit_pf');            
    makedir(fit_dir);

    for n=1:length(data)
        fname_subj = fullfile(fit_dir, sprintf('subj_%04d.mat', n));
        if ~exist(fname_subj, 'file')
            dat = struct('outcome', data{n}.bird, 'action', data{n}.bucket);
            
            fit_model(fname_subj, config, dat, model, opt_var, InitialX);
        end
        f = load(fname_subj);    
        loglik(n) = f.pf_observed.loglik;
        parameters(n,:) = f.pf_observed.x;
        dynamics_mean = f.pf_observed.signal.dynamics_mean;
        dynamics{n} = struct('val', dynamics_mean.val, 'vol', dynamics_mean.vol, 'sto', dynamics_mean.sto, 'lr', dynamics_mean.lr);
    end

    logf = -.5*loglik -0.5*log(2*pi);
    num_observations = numel(data{1}.bucket);
    num_parameters = width(InitialX);
    lme = logf -0.5*num_parameters*log(num_observations);
    
    save(fname, 'loglik', 'parameters', 'dynamics', 'lme');    
end

end

% -------------------------------------------------------------------------
function bayesopt_results = fit_model(fname, config, data, model, opt_var, InitialX)


fun = @(x)model(x, data.outcome, 'fit');
opt_fun = @(x)optimization(x, fun, config.num_inits, data.action, 'rmean');

rng(config.seed);
bayesopt_results = bayesopt(opt_fun, opt_var,'IsObjectiveDeterministic', false, 'Verbose', 1,...
    'PlotFcn',{},'AcquisitionFunctionName','expected-improvement-plus',...
    'InitialX',InitialX, ...
    'MaxObjectiveEvaluations', config.num_bayespoint,'NumSeedPoints', config.num_seeds);


[pf_observed, pf_estimated] = fit_post(bayesopt_results, data, model, config);

if ~isempty(fname)
    save(fname, 'bayesopt_results', 'config', 'pf_estimated', 'pf_observed');
end

end


function [pf_observed, pf_estimated] = fit_post(bayesopt_results, dat, model, config)
        
fun = @(x)model(x, dat.outcome, 'fit');
opt_fun = @(x)optimization(x, fun, config.num_inits, dat.action, 'mean');

rng(config.seed);
[loglik, ~, log_likelihood] = opt_fun(bayesopt_results.XAtMinObjective);
pf_observed = struct('x', bayesopt_results.XAtMinObjective, 'loglik', loglik,...
                'log_likelihood',log_likelihood);    
rng(config.seed);
pf_observed.signal =  make_signal(pf_observed.x, model, dat, config.num_inits);


rng(config.seed);
[loglik, ~, log_likelihood] = opt_fun(bayesopt_results.XAtMinEstimatedObjective);
pf_estimated = struct('x', bayesopt_results.XAtMinEstimatedObjective, 'loglik', loglik,...
             'log_likelihood',log_likelihood);    
rng(config.seed);
pf_estimated.signal =  make_signal(pf_estimated.x, model, dat, config.num_inits);


end

function [objective, c, loglik, val] = optimization(x, fun, num_init, action, obj_stats)
ll = zeros(1,num_init);
log_likelihood = zeros(num_init, 4);
val = cell(1,num_init);
for k = 1:num_init
    val{k} = fun(x);
    log_likelihood(k,:) = sum( ((val{k}-action).^2));    
    ll(k) = sum(log_likelihood(k,:));
end
if strcmp(obj_stats, 'mean')
    objective = mean(ll);
elseif strcmp(obj_stats, 'rmean')
    objective = mean(ll)+std(ll)/sqrt(num_init);
else
    error('unknown obj_stats!');
end

c = [];
loglik = mean(ll);
end

function sig = make_signal(parameters, model, dat, num_inits)

for i=1:num_inits
    [~, vars_i] = model(parameters, dat.outcome, 'sim');
    snames = fieldnames(vars_i);
    for j = 1:length(snames)
        if i==1
            dynamics_mean.(snames{j}) = 0;
        end
        dynamics_mean.(snames{j}) = dynamics_mean.(snames{j}) + vars_i.(snames{j})/num_inits;
    end  
end
sig = struct('dynamics_mean', dynamics_mean, 'parameters', parameters);
end
