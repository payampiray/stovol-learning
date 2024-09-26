function [val, vars] = model_weber(parameters, observations, mode)

if strcmp(mode, 'config')
    pnames = {'alpha', 'zeta', 'sigma'};    
    mu1_var = optimizableVariable(pnames{1}, [0.01, 1]);
    mu2_var = optimizableVariable(pnames{2}, [0.01, 1]);
    sigma_var = optimizableVariable(pnames{3}, [0.05, 1]);
    opt_var = [mu1_var, mu2_var, sigma_var];
    InitialX = array2table([0.01, 0.01, .05; 1, 1, .05; 1, 1, 1], 'VariableNames', pnames);  
    InitObj = [];

    val = nan;
    vars = struct('opt_var', opt_var, 'InitialX', InitialX, 'InitObj', InitObj, 'parameters_name', {pnames});
    return;
end

[val, vars] = pf_process(parameters, observations);

end

function [q] = transition_func(q, delta, params)
alpha = params.alpha;
zeta = params.zeta;

epsil = randn(size(q)).*abs(delta)*zeta;

q = q + alpha*delta + epsil;

end

%--------------------------------------------------------------------------
function [val, vars] = pf_process(parameters, observations, config)
if nargin < 3
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
    [val(:,i)] = pf_process_block(parameters, observations(:,i), config);
    vars = [];
end

end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [val, vars] = pf_process_block(parameters, observations, config)
num_particles = config.num_particles;

weights = 1/num_particles*ones(1, num_particles);

N = length(observations);

m = 60 + zeros(1, num_particles);

val_post = nan(N, 1);
val = nan(N, 1);
lr = nan(N, 1);
for t=1:size(observations,1)    
    val(t) = sum(m.*weights);    
    
    delta = observations(t) - m;
    m_post = transition_func(m, delta, parameters);

    a = (m_post - m)./delta;
    m = m_post;
    
    likelihood = normpdf(observations(t), m, parameters.sigma);
    [idx, weights] = resample(likelihood, weights, config);    
    m = m(idx);
    
    lr(t) = sum(a.*weights);           
    val_post(t) = sum(m.*weights);        
end

vars = struct('val',val,'lr', lr, 'val_post', val_post);
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

end

function [idx, w] = resampling(resampling_strategy, w)
N  = length(w);
switch resampling_strategy
    case 'systematic'
        edges = min([0 cumsum(w)],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        u = rand/N;
        [~, idx] = histc(u:1/N:1, edges); %#ok<HISTC> 
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

