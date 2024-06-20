function [tbl] = other_rw(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),'model_rw.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    N = length(data);
    config = struct('bound', [0 0 0 0; 1 1 1 1]);
    loglik = nan(N, 1);
    lme = nan(N, 1);
    for n=1:N        
        [parameters(n, :), loglik(n), lme(n)] = tools_fit(data{n}, @model_RW, config);
        fprintf('%03d\n', n);
    end
        
    save(fname, 'parameters', 'loglik', 'lme');      
end
f = load(fname);

parameters = f.parameters;
x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'alpha1','alpha2','alpha3','alpha4'};
tbl.columns = {'25%', '50%', '75%'};

end

% -------------------------------------------------------------------------
function [objective] = model_RW(f, outcome, action)
alpha = f';

ndim = size(outcome,2);
N = size(outcome,1);
val = nan(N,ndim);

m = 60+zeros(1,ndim);
for t=1:N    
    val(t, :) = m;
    m = m + alpha.*(outcome(t, :)-m);
end

log_likelihood = sum(-0.5*(val-action).^2 -0.5*log(2*pi));
loglik = sum(log_likelihood);

objective = -loglik;
end