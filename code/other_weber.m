function [tbl] = other_weber(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment), 'model_weber.mat');
    
if ~exist(fname, 'file')
    data = get_data(experiment);
    config = struct('num_bayespoint', 200, 'num_seeds', 5, 'seed', 0, 'num_inits', 10);            
    model = @model_weber;    
    
    [~, vars] = model([], [], 'config');
    opt_var = vars.opt_var;
    InitialX = vars.InitialX;    

    loglik = nan(length(data), 1);
    dynamics = cell(length(data), 1);

    fit_dir = fullfile('..',sprintf('experiment%d', experiment), 'fit_weber');            
    mkdir(fit_dir);

    for n=1:length(data)
        fname_subj = fullfile(fit_dir, sprintf('subj_%04d.mat', n));        

        if ~exist(fname_subj, 'file')
            dat = struct('outcome', data{n}.bird, 'action', data{n}.bucket);
            
            fit2_fit(fname_subj, config, dat, model, opt_var, InitialX);
        end
        f = load(fname_subj);    
        loglik(n) = f.pf_estimated.loglik;
        parameters(n,:) = f.pf_observed.x;
        dynamics_mean = f.pf_estimated.signal.dynamics_mean;
        dynamics{n} = struct('val', dynamics_mean.val, 'lr', dynamics_mean.lr,  'val_post', dynamics_mean.val_post);
    end

    logf = -.5*loglik -0.5*log(2*pi);
    num_observations = numel(data{1}.bucket);
    num_parameters = width(InitialX);
    lme = logf -0.5*num_parameters*log(num_observations);
    
    save(fname, 'loglik', 'parameters', 'dynamics', 'lme');

end

f = load(fname);
parameters = f.parameters.Variables;

x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'alpha', 'zeta', 'sigma'};
tbl.columns = {'25%', '50%', '75%'};
end
