function kf_recovery

[fit_effect, true_effect] = do_run();

c = corr(fit_effect, true_effect);
c = diag(c);

me = median(abs(fit_effect - true_effect));
se = se_median(abs(fit_effect - true_effect));

figure;

for i=1:4
subplot(2,2,i);
scatterplot(fit_effect(:,i), true_effect(:, i))
    
end
end

function [fit_effect, true_effect] = do_run()
fname = fullfile('..', 'analysis', 'fit_kf_model1', 'sim_4params.mat');

nsim = 100;
noise_action = 0;
num_parameters = 4;
init_range = 50;

config.init_range = init_range;
config.num_parameters = num_parameters;

data = get_data(1);

A = ([0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1;1 -1 -1 1])^-1;
A_inv = [0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1;1 -1 -1 1];

fit_effect = nan(nsim, num_parameters);
fit_kalman_parameter = nan(nsim, num_parameters);
true_effect = nan(nsim, num_parameters);
true_kalman_parameter = nan(nsim, num_parameters);


rng(0);
true_lambda = rand(num_parameters, nsim)*init_range;


if ~exist(fname, 'file')
    for n=1:nsim
        outcome = data{n}.bag;

        lambda = true_lambda(:, n);
        true_parameters = A_inv*lambda;
        
        [action] = model_kalman(true_parameters, A, outcome, noise_action);
        sim = struct('bag', outcome, 'bucket', action);
        fit(n) = fit1_fit(sim, config); %#ok<AGROW> 
    
        fit_effect(n, :) = fit(n).effect;
        fit_kalman_parameter(n, :) = fit(n).kalman_parameter;
        true_effect(n, :) = true_parameters;
        true_kalman_parameter(n, :) = lambda;    
    end
    save(fname, 'fit', 'fit_effect', 'fit_kalman_parameter', 'true_effect', 'true_kalman_parameter');    
end
f = load(fname);
fit_effect = f.fit_effect;
true_effect = f.true_effect;
end

function [action] = model_kalman(f, M, outcome, noise_action)
lambda = (M*f)';

ndim = size(outcome,2);
N = size(outcome,1);
val = nan(N,ndim);

m = 60+zeros(1,ndim);
r = 1*ones(1,ndim);

for t=1:N    
    val(t, :) = m;
    a = (r+lambda)./(r+1+lambda);
    m = m + a.*(outcome(t, :)-m);
    r = (1-a).*(r+lambda);        
end

action = val + noise_action*randn(size(val));

end