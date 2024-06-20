function tbl = other_hgf(experiment)
if nargin<1, experiment = 1; end
fname = fullfile('..',sprintf('experiment%d', experiment), 'model_hgf.mat');

data = get_data(experiment);

if ~exist(fname, 'file')
    N = length(data);
    lb = [0 0 1];
    ub = [4 .5 100];
    config = struct('bound', [lb; ub]);
    loglik = nan(N, 1);
    lme = nan(N, 1);    
    dynamics = cell(N, 1);    
    for n=1:N        
        [parameters(n, :), loglik(n), lme(n)] = tools_fit(data{n}, @model_hgf, config);
        [~, dynamics{n}] = model_hgf(parameters(n, :), data{n}.bag, data{n}.bucket);
        fprintf('%03d\n', n);
    end
        
    save(fname, 'parameters', 'loglik', 'lme', 'dynamics');    
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
