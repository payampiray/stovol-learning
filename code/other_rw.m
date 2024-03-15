function [tbl, parameter] = other_rw(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),'model_rw.mat');

data = get_data(experiment);
config.num_parameters = 4;

if ~exist(fname, 'file')
    num_params = config.num_parameters;
    N = length(data);
    lme = nan(N, 1);
    parameter = nan(N, num_params);
    for n=1:N        
        fit(n) = fit1_rl(data{n}, config); %#ok<AGROW> 
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
function [fit] = fit1_rl(dat, config)
if isempty(config)
    config = struct('num_parameters', 4);
end

p = inputParser;
p.addParameter('num_parameters', 4);
p.addParameter('init_range', 1);
p.addParameter('bound', 1);
p.addParameter('num_init', 10);

p.parse(config);
config    = p.Results;


action = dat.bucket;
outcome = dat.bag;


options = optimoptions('fmincon','Display','off');
x0 = .5*[1 1 1 1]';
b = zeros(4,1);
A = eye(4,4);

lb = [0 0 0 0];
ub = [1 1 1 1];
l0 = rand(4, config.num_init-1)*config.init_range;
x0 = [x0 l0];

lb = -config.bound*lb';
ub = +config.bound*ub';
init_range = config.init_range;
num_init = config.num_init;

l0 = rand(4, config.num_init-1)*config.init_range;
x0 = [x0 l0];

fun = @(f)model_rl(f, outcome, action);

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

config = struct('lb', lb, 'ub', ub, 'init_range', init_range, 'num_init', num_init, 'x0', x0);
optim = struct('flag', flag, 'gradient', g, 'hessian', H, 'is_H_pos', is_H_pos, 'neg_loglik_vec', neg_loglik_vec, 'options', options);

fit = struct('lme', lme, 'parameter', x', 'loglik', loglik, 'Ainvdiag', {{Ainvdiag}}, 'config', config, 'optim', optim);
end

function [objective] = model_rl(f, outcome, action)
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