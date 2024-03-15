function tbl = other_hgf(experiment)
if nargin<1, experiment = 1; end
fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hgf.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    num_params = 3;
    N = length(data);
    lme = nan(N, 1);
    parameters = nan(N, num_params);
    dynamics = cell(N, 1);
    for n=1:N        
        fit(n) = fit_hgf(data{n}); %#ok<AGROW> 
        lme(n) = fit(n).lme;
        [~, dynamics{n}] = model_hgf(fit(n).parameter, data{n}.bag, data{n}.bucket);
        parameters(n, :) = fit(n).parameter;
        fprintf('%03d\n', n);
    end
        
    save(fname, 'fit', 'lme', 'parameters', 'dynamics');    
end
f = load(fname);
parameters = f.parameters;
parameters(:, 2) = log(parameters(:, 2));
x = prctile(parameters, [25 50 75])';
tbl.data = x;
tbl.rows = {'nu','log(kappa)','hgf alpha'};
tbl.columns = {'25%', '50%', '75%'};

end

% -------------------------------------------------------------------------
function [fit] = fit_hgf(dat)

config = struct('init_range', 1, 'num_init', 10);

action = dat.bucket;
outcome = dat.bag;


options = optimoptions('fmincon','Display','off');

x0 = [2 .25 50]';
b = zeros(3,1);
A = eye(3);

lb = [0 0 1]';
ub = [4 .5 100]';
num_init = config.num_init;

l0 = lb + (ub - lb).*rand(3, config.num_init-1);
x0 = [x0 l0];

fun = @(f)model_hgf(f, outcome, action);

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
parameter = x';

fit = struct('lme', lme, 'parameter', parameter, 'loglik', loglik, 'Ainvdiag', {{Ainvdiag}}, 'config', config, 'optim', optim);
end

function [objective, signals] = model_hgf(f, y, action)
nu = f(1);
kappa = f(2);
expomega = 1;
alpha = f(3);

[N,Q]  = size(y);
y      = [zeros(1,Q); y];
mu3    = nan(N+1,Q);
mu2    = nan(N+1,Q);
sigma2 = nan(N+1,Q);
sigma3 = nan(N+1,Q);

mu2(1,:)    = (log(alpha))/kappa;
sigma2(1,:) = 1;
mu3(1,:)    = 1;
sigma3(1,:) = 1;

LR     = nan(N+1,Q);
v      = nan(N+1,Q);
v(1,:) = expomega*exp(kappa*mu3(1,:));

bad_sigma = 0;
for n  = 2:(N+1)

    expmu3          = expomega*exp(kappa*mu3(n-1,:));
    
    delta1          = y(n,:)- mu2(n-1,:);                                  
    LR(n,:)         = (expmu3+sigma2(n-1,:))./(expmu3+sigma2(n-1,:) + alpha);
    mu2(n,:)        = mu2(n-1,:) + LR(n,:).*delta1;                         % Eq 51
    sigma2(n,:)     = LR(n,:)*alpha;                                        % Eq 50
        
    pihat3          = (sigma3(n-1,:) + nu).^-1;                             % Eq 31
    w2              = expmu3 ./ (expmu3 + sigma2(n-1,:) );                  % Eq 32
    r2              = (expmu3 - sigma2(n-1,:))./(expmu3 + sigma2(n-1,:) );  % Eq 33
    delta2          = (sigma2(n,:) +...
        (mu2(n,:)-mu2(n-1,:)).^2)./(sigma2(n-1,:) + expmu3) - 1;            % Eq 34
    
    pi3             = pihat3 + (kappa^2/2).*w2.*(w2+r2.*delta2);            % Eq 29    
    if any(pi3 <= 0)
        bad_sigma = true;
    end
    sigma3(n,:)     = pi3.^-1;    
    mu3(n,:)        = mu3(n-1,:) + (kappa/2).*sigma3(n,:).*w2.*delta2;      % Eq 30
      
    
    v(n,:)     = expomega*exp(kappa*mu3(n,:));
end

mu2(1, :) = [];
v(1, :) = [];
LR(1, :) = [];
signals = struct('val',mu2,'volatility',v,'lr',LR);

log_likelihood = sum(-0.5*(mu2-action).^2 -0.5*log(2*pi));
loglik = sum(log_likelihood);

objective = -loglik;

if bad_sigma
    objective = 10^6;
end

end
