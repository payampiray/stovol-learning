function [tbl] = other_dbd4(experiment)
if nargin<1, experiment = 1; end
fname = fullfile('..',sprintf('experiment%d', experiment), 'model_dbd4.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    N = length(data);
    lme = nan(N, 1);
    parameters = nan(N, 4);
    dynamics = cell(N, 1);
    for n=1:N        
        fit(n) = fit_model(data{n}); %#ok<AGROW> 
        lme(n) = fit(n).lme;
        parameters(n, :) = fit(n).parameters;
        [~, dynamics{n}] = model_dbd(fit(n).parameters', data{n}.bag, data{n}.bucket);
        fprintf('%03d\n', n);
    end
        
    save(fname, 'fit', 'lme', 'parameters', 'dynamics');    
end
f = load(fname);
parameters = f.parameters;
x = prctile(log(parameters), [25 50 75])';
tbl.data = x;
tbl.rows = {'theta1','theta2','theta3','theta4'};
tbl.columns = {'25%', '50%', '75%'};

end

% -------------------------------------------------------------------------
function [fit] = fit_model(dat)
num_init = 10;
init_range = 1;

action = dat.bucket;
outcome = dat.bag;


options = optimoptions('fmincon','Display','off');

x0 = .5*[1 1 1 1]';
b = zeros(4,1);
A = eye(4,4);

lb = [0 0 0 0];
ub = [1 1 1 1];
l0 = rand(4, num_init-1)*init_range;
x0 = [x0 l0];


l0 = rand(4, num_init-1);
x0 = [x0 l0];

fun = @(f)model_dbd(f, outcome, action);

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

fit = struct('lme', lme, 'parameters', x', 'loglik', loglik, 'Ainvdiag', {{Ainvdiag}});
end


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

