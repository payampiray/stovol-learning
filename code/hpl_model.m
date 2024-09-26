function [val, vars] = hpl_model(parameters, observations, mode)

if strcmp(mode, 'config')
    pnames = {'mu_s', 'mu_v', 'sigma'};    
    mu_s = optimizableVariable(pnames{1}, [0.01, 1]);
    mu_v = optimizableVariable(pnames{2}, [0.01, 1]);
    sigma_var = optimizableVariable(pnames{3}, [.05, 1]);
    opt_var = [mu_s, mu_v, sigma_var];
    InitialX = array2table([0.01, 0.01, .05; 1, 1, .05;1, 1, 1], 'VariableNames', pnames);  
    rand_func = @(n)[.01+.99*rand(n,1), .01+.99*rand(n,1), .05+.95*rand(n,1)];

    val = nan;
    vars = struct('opt_var', opt_var, 'InitialX', InitialX, 'parameters_name', {pnames}, 'rand_func', rand_func);
    return;
end

[val, vars] = pf_process(parameters, observations, mode);

end

function [v, s, mix] = transition_func(v, s, params)
mu_s = params.mu_s;
mu_v = params.mu_v;
sigma = params.sigma;

mix_s = binornd(1, mu_s, size(s));
mix_v = binornd(1, mu_v, size(v));
mix = [mix_s; mix_v];

e1 = sigma*randn(size(s));
e2 = sigma*randn(size(s));

s = mix_s.*(log(s)+e1) + (1-mix_s).*log(s);
v = mix_v.*(log(v)+e2) + (1-mix_v).*log(v);

s = exp(s);
v = exp(v);
end

%--------------------------------------------------------------------------
function [val, vars] = pf_process(parameters, observations, mode, config)
if nargin < 4
    config.resampling_strategy = 'systematic';
end

p = inputParser;
p.addParameter('resampling_strategy', 'systematic');
p.addParameter('resample_percentaje', .5);
p.addParameter('num_particles', 10000);
p.parse(config);
config    = p.Results;

val = nan(size(observations)); 
for i=1:size(observations,2)
    [val(:,i), vars_i] = pf_process_block(parameters, observations(:,i), config);
    if strcmp(mode, 'fit')        
        vars = [];
    elseif strcmp(mode, 'sim')
        snames = fieldnames(vars_i);
        for j = 1:length(snames)
            if i==1
                vars.(snames{j}) = vars_i.(snames{j});
            else
                vars.(snames{j})(:,i) = vars_i.(snames{j});
            end
        end
    end
end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [val, vars] = pf_process_block(parameters, observations, config)
num_particles = config.num_particles;

weights = 1/num_particles*ones(1, num_particles);

N = length(observations);

u = 1;
w = 10;
s = u*ones(1,num_particles);
v = u*ones(1,num_particles);
w = w*ones(1,num_particles);

m = 60 + zeros(1, num_particles);
w = w.*ones(1, num_particles);

val = nan(N, 1);
vol = nan(N, 1);
sto = nan(N, 1);
lr = nan(N, 1);

for t=1:size(observations,1)    
    val(t) = sum(m.*weights);    
    
    [v, s] = transition_func(v, s, parameters);

    vol(t) = sum(v.*weights);
    sto(t) = sum(s.*weights);    
    
    likelihood = normpdf(observations(t),m, sqrt(v+s+w));
    [idx, weights] = resample(likelihood, weights, config);    
    s = s(idx);
    v = v(idx);
    m = m(idx);
    w = w(idx);
        
    [m, w, k] = kalman(observations(t),m, w, s, v);
    lr(t) = sum(k.*weights);
end

vars = struct('val',val,'lr', lr, 'vol', vol, 'sto', sto);
end

function [idx, weights, resampled] = resample(likelihood, weights, config)
NumParticles = length(likelihood);

weights = weights.*(likelihood+eps);
weights = weights/sum(weights);

Neff = 1/sum(weights.^2);
resample_percentaje = config.resample_percentaje;
Nt = resample_percentaje*NumParticles;
idx = 1:NumParticles;
resampled = 0;
if Neff < Nt
    [idx, weights] = resampling(config.resampling_strategy, weights);
    resampled = 1;
end

if any(isnan(weights))
    disp('')
end

end

function [idx, w] = resampling(resampling_strategy, w)
N  = length(w);
switch resampling_strategy
    case 'systematic'
        edges = min([0 cumsum(w)],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        u = rand/N;
        [~, idx] = histc(u:1/N:1, edges);
    case 'multinomial'        
        u = rand(N,1);
        wc = cumsum(w');
        wc = wc(:);
        wc = wc/wc(N);
        [~,ind1] = sort([u;wc]);
        ind2 = find(ind1<=N);
        idx = ind2'-(0:N-1);
   otherwise
      error('Resampling strategy is unknown');
end

w = ones(size(w))/N;

end

function [m, w, k, delta]=kalman(outcome, m, w, s, v)
delta = outcome - m;
k = (w+v)./(w+v+s);
m = m + k.*(outcome-m);
w = (1-k).*(w+v);
end
