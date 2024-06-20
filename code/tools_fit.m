function [parameters, loglik, lme] = tools_fit(dat, model, config)

num_parameters = size(config.bound, 2);
num_init = 10;

action = dat.bucket;
outcome = dat.bag;


x0 = mean(config.bound, 1)';
b = zeros(num_parameters, 1);
A = eye(num_parameters);


lb = config.bound(1, :);
ub = config.bound(2, :);
l0 = rand(num_parameters, num_init-1).*(ub - lb)';
x0 = [x0 l0];

fun = @(f)model(f, outcome, action);

x = cell(1, num_init);
neg_loglik_vec = nan(1, num_init);
flag = cell(1, num_init);
g = cell(1, num_init);
H = cell(1, num_init);
options = optimoptions('fmincon','Display','off');

for i=1:num_init
    [x{i}, neg_loglik_vec(i), flag{i}, ~,~,g{i}, H{i}] = fmincon(fun, x0(:, i), -A, b,[],[], lb, ub,[],options);
end

[~, i] = min(neg_loglik_vec);
neg_loglik = neg_loglik_vec(i);
x = x{i};
H = H{i};

[~,is_H_pos] = chol(H);
is_H_pos = ~logical(is_H_pos);
if ~is_H_pos
    x = x0(:, 1);
    neg_loglik = fun(x);
end

loglik = -neg_loglik;
parameters = x';
lme = loglik -.5*num_parameters*log(200);
end