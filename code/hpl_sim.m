function hpl_sim
experiment = 1;
% num_sim = 100;

nn = 1:200;

fname = fullfile('..',sprintf('experiment%d', experiment), 'hpl_sim.mat');   

if ~exist(fname, 'file')
    data = get_data(experiment);
    config = struct('num_bayespoint', 200, 'num_seeds', 5, 'seed', 0, 'num_inits', 10);    
    model = @model_pf;    
    
    [~, vars] = model([], [], 'config');
    opt_var = vars.opt_var;
    InitialX = vars.InitialX;

    sim_data = cell(length(data), 1);
    fit_dir = fullfile(sprintf('experiment%d', experiment), 'hpl_sims');
    mkdir(fit_dir);

    for n=nn
        fname_subj = fullfile(fit_dir, sprintf('subj_%04d.mat', n));

        if ~exist(fname_subj, 'file')
            [dat, x_true] = simulate(model, config, n, data{n}.bag);            
            fit2_fit(fname_subj, config, dat, model, opt_var, InitialX);
            f = load(fname_subj);
            f.true_parameters = x_true;
            f.sim = dat;            
            save(fname_subj, '-struct', 'f')
        end
        f = load(fname_subj);    
        sim = f.sim;
        sim_data{n} = struct('outcome', sim.outcome, 'action', sim.action);
        fitted_parameters(n,:) = f.pf_estimated.x;
        true_parameters(n,:) = f.true_parameters;
    end
    save(fname, 'fitted_parameters', 'true_parameters', 'sim_data');
end

end

function [sim, true_parameters] = simulate(model, config, n, outcome)
[~, vars] = model([], [], 'config');
pnames = vars.parameters_name;
rand_func = vars.rand_func;

rng(n);
true_parameters = rand_func(1);

true_parameters = array2table(true_parameters, 'VariableNames', pnames);

num_inits = config.num_inits;
action = 0;
for k = 1:num_inits
    [val] = model(true_parameters, outcome, 'sim');    
    action = action + val/num_inits;
end

sim = struct('outcome', outcome, 'action', action);

end