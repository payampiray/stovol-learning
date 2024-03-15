function [fx, labels, kalman_parameter] = kf_fit(experiment)
if nargin<1, experiment = 1; end

fname = fullfile('..',sprintf('experiment%d', experiment),'model_kf4.mat');

if ~exist(fname, 'file')
    data = get_data(experiment);    
    num_params = 4;
    N = length(data);
    lme = nan(N, 1);
    effect = nan(N, num_params);
    kalman_parameter = nan(N, num_params);
    dynamics = cell(N, 1);
    for n=1:N        
        kalman(n) = fit_model(data{n}); %#ok<AGROW> 
        lme(n) = kalman(n).lme;
        effect(n, :) = kalman(n).effect;
        kalman_parameter(n, :) = kalman(n).kalman_parameter;
        dynamics{n} = kalman_model(kalman(n).kalman_parameter, data{n}.bag);
        fprintf('%03d\n', n);
    end

    save(fname, 'kalman', 'lme', 'effect', 'kalman_parameter', 'dynamics');
end

f = load(fname);
fx = f.effect;
fx = fx(:, [2:4 1]);
labels = {'lambda_s', 'lambda_v', 'lambda_i', 'lambda_m'};
kalman_parameter = f.kalman_parameter;
end

function fit = fit_model(dat)
config = struct('num_parameters', 4, 'init_range', 50, 'bound', 100, 'num_init', 10);

action = dat.bucket;
outcome = dat.bag;


options = optimoptions('fmincon','Display','off');

if config.num_parameters == 4
    x0 = [1 0 0 0]';
    b = zeros(4,1);
    A = ([0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1;1 -1 -1 1])^-1;
    A_inv = [0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1;1 -1 -1 1];

    lb = [0 1 1 1];
    ub = [1 1 1 1];
    l0 = rand(4, config.num_init-1)*config.init_range;
    x0 = [x0 (A^-1)*l0];

elseif config.num_parameters == 3
    x0 = [1 0 0]';
    b = zeros(4,1);
    A = ([0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1;1 -1 -1 1])^-1;
    A(:, end) = [];
    A_inv = [0.25 0.25 0.25 0.25;1 1 -1 -1;-1 1 -1 1];

    lb = [0 1 1];
    ub = [1 1 1];

elseif config.num_parameters == 1
    x0 = 1;
    b = 0;
    A = 1;
    A_inv = [0.25 0.25 0.25 0.25];

    lb = 0;
    ub = 1;    
end

lb = -config.bound*lb';
ub = +config.bound*ub';
init_range = config.init_range;
num_init = config.num_init;

l0 = rand(4, config.num_init-1)*config.init_range;
x0 = [x0 A_inv*l0];

fun = @(f)model_kalman(f, A, outcome, action);


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

config = struct('A',A,'b',b,'lb', lb, 'ub', ub, 'init_range', init_range, 'num_init', num_init, 'x0', x0);
optim = struct('flag', flag, 'gradient', g, 'hessian', H, 'is_H_pos', is_H_pos, 'neg_loglik_vec', neg_loglik_vec, 'options', options);
lambda = (A*x)';

fit = struct('lme', lme, 'effect', x', 'kalman_parameter', lambda, 'loglik', loglik, 'Ainvdiag', {{Ainvdiag}}, 'config', config, 'optim', optim);
end


function [objective] = model_kalman(f, M, outcome, action)
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

log_likelihood = sum(-0.5*(val-action).^2 -0.5*log(2*pi));
loglik = sum(log_likelihood);

objective = -loglik;
end
