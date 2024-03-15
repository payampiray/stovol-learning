function [tbl, parameter] = other_kfem(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),'model_kfem.mat');

data = get_data(experiment);
config.num_parameters = 2;

if ~exist(fname, 'file')
    num_params = config.num_parameters;
    N = length(data);
    lme = nan(N, 1);
    parameter = nan(N, num_params);
    for n=1:N        
        fit(n) = fit1_kfem(data{n}, config); %#ok<AGROW> 
        lme(n) = fit(n).lme;
        parameter(n, :) = fit(n).parameter;
        fprintf('%03d\n', n);
    end
        
    save(fname, 'fit', 'lme', 'parameter');    
end
f = load(fname);
lme = f.lme;
parameters = f.parameter;
x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'alpha1','alpha2','alpha3','alpha4'};
tbl.columns = {'25%', '50%', '75%'};

end

% -------------------------------------------------------------------------
function [fit] = fit1_kfem(dat, config)
if isempty(config)
    config = struct('num_parameters', 2);
end

p = inputParser;
p.addParameter('num_parameters', 2);
p.addParameter('init_range', 1);
p.addParameter('bound', 1);
p.addParameter('num_init', 10);

p.parse(config);
config    = p.Results;


action = dat.bucket;
outcome = dat.bag;


options = optimoptions('fmincon','Display','off');
x0 = 50*[1 1]';
b = zeros(2,1);
A = eye(2, 2);

lb = [0 0];
ub = [1 1]*100;
l0 = rand(2, config.num_init-1)*config.init_range;
x0 = [x0 l0];

num_init = config.num_init;

fun = @(f)model_kfem(f, outcome, action);

x = cell(1, num_init);
neg_loglik_vec = nan(1, num_init);
flag = cell(1, num_init);
g = cell(1, num_init);
H = cell(1, num_init);

for i=1:num_init
    [x{i}, neg_loglik_vec(i), flag{i}, ~,~,g{i}, H{i}] = fmincon(fun, x0(:, i), -A, b,[],[], lb, ub,[],options);
end

[~, i] = min(neg_loglik_vec);
neg_loglik = neg_loglik_vec(i);
x = x{i};
flag = flag{i};
g = g{i};
H = H{i};

[~,is_H_pos] = chol(H);
is_H_pos = ~logical(is_H_pos);
if ~is_H_pos
    x = x0(:, 1);
    neg_loglik = fun(x);
%     H = prior.precision;    
end

loglik = -neg_loglik;

num_params = length(x0);
Ainvdiag = diag(H);
lme = loglik -.5*num_params*log(200);

config = struct('lb', lb, 'ub', ub, 'num_init', num_init, 'x0', x0);
optim = struct('flag', flag, 'gradient', g, 'hessian', H, 'is_H_pos', is_H_pos, 'neg_loglik_vec', neg_loglik_vec, 'options', options);

fit = struct('lme', lme, 'parameter', x', 'loglik', loglik, 'Ainvdiag', {{Ainvdiag}}, 'config', config, 'optim', optim);
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
