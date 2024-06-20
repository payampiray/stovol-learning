function tbl = other_kfem(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),'model_kfem.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    N = length(data);
    config = struct('bound', [0 0; 100 100]);
    loglik = nan(N, 1);
    lme = nan(N, 1);    
    for n=1:N        
        [parameters(n, :), loglik(n), lme(n)] = tools_fit(data{n}, @model_kfem, config);
        fprintf('%03d\n', n);
    end
        
    save(fname, 'parameters', 'loglik', 'lme');    
end
f = load(fname);

parameters = f.parameters;
x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'beta_v','beta_s'};
tbl.columns = {'25%', '50%', '75%'};

end

function [objective] = model_kfem(f, outcome, action)
alpha_v = 1;
alpha_s = 1;
beta_v = f(1);
beta_s = f(2);

ndim = size(outcome,2);
N = size(outcome,1);
val = nan(N,ndim);

m = 60+zeros(1,ndim);
w = 10*ones(1,ndim);
v = 0*ones(1,ndim);
s = 0*ones(1,ndim);
m_pre = m;
w_pre = w;
a = 1;

for t=1:N
    
    val(t, :) = m;           

    e = (m-m_pre).^2 + w + w_pre -2*w_pre.*(1-a) + beta_v;

    lambda = 1/(t-1 + alpha_v);
    v = v + lambda*(e -v);

    e = (outcome(t,:)-m).^2 + w + beta_s;
    lambda = 1/(t + alpha_s);
    s = s + lambda*(e -s);

    m_pre = m;
    w_pre = w;
    a = (w+v)./(w+s+v);
    m = m + a.*(outcome(t, :)-m);
    w = (1-a).*(w+v);    

    lr(t, :) = a;
end

log_likelihood = sum(-0.5*(val-action).^2 -0.5*log(2*pi));
loglik = sum(log_likelihood);

objective = -loglik;

end
