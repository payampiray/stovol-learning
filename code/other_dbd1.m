function [tbl] = other_dbd1(experiment)
if nargin<1, experiment = 1; end
fname = fullfile('..',sprintf('experiment%d', experiment), 'model_dbd1.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    N = length(data);
    config = struct('bound', [0; 1]);
    parameters = nan(N, 1);
    loglik = nan(N, 1);
    lme = nan(N, 1);    
    dynamics = cell(N, 1);    
    for n=1:N        
        [parameters(n, :), loglik(n), lme(n)] = tools_fit(data{n}, @model_dbd, config);
        [~, dynamics{n}] = model_dbd(parameters(n, :), data{n}.bag, data{n}.bucket);
        fprintf('%03d\n', n);
    end
        
    save(fname, 'parameters', 'loglik', 'lme', 'dynamics');      
end
f = load(fname);

parameters = f.parameters;
x = prctile((parameters), [25 50 75]);
tbl.data = x;
tbl.rows = {'theta'};
tbl.columns = {'25%', '50%', '75%'};

end

% -------------------------------------------------------------------------
function [objective, dynamics] = model_dbd(x, outcome, action)
theta = x';

ndim = size(outcome,2);
N = size(outcome,1);
val = nan(N,ndim);

beta = zeros(1,4);
h = ones(1,4);

m = 60+zeros(1,ndim);
for t=1:N    
    val(t, :) = m;
    delta = outcome(t, :) - m;
    m_old = m;

    beta = beta + theta.*delta.*h;
    alpha = exp(beta);
    m = m + alpha.*delta;
    h = h + subplus(1-alpha) + alpha.*delta;

    beta_val(t, :) = beta;
    lr(t, :) = alpha;
end

log_likelihood = sum(-0.5*(val-action).^2 -0.5*log(2*pi));
loglik = sum(log_likelihood);

objective = -loglik;

dynamics = struct('val', val, 'lr', lr, 'beta', beta_val);
end

