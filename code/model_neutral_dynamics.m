function model_neutral_dynamics(data_set, nr, nc, subplots)
if nargin<1
    data_set = 1;
end

mode = 2;
wind = 5;
nsim = 100;
makedir(fullfile('..','glm'));
fname = fullfile('..','glm',sprintf('generate_rev_data%d_nsim%d_m%dw%d.mat', data_set, nsim, mode, wind));

[data] = get_data(data_set);

if ~exist(fname, 'file')
    [xb, xe] = rev_glm_get_task_index(data, mode);

    for n=1:nsim
        rng(n);
        data_pf{n} = sim_model(data_set);
        [pb(:,:,n), pe(:,:,n)] = rev_glm_get_task_index(data_pf{n}, mode);
        
        rng(n);
        data_kf{n} = sim_kalman(data_set);
        [kb(:,:,n), ke(:,:,n)] = rev_glm_get_task_index(data_kf{n}, mode);
    
        fprintf('%04d is done\n', n)
    end
    save(fname,'data_pf', 'data_kf', 'pb', 'pe', 'kb', 'ke', 'xb', 'xe');
end

f = load(fname);
xb = f.xb;
xe = f.xe;
pb = mean(f.pb, 3);
kb = mean(f.kb, 3);
% pe = serr(f.pb, 3);



cp = corr(xb, pb, 'type', 'spearman');
cp0 = corr(xb(:), pb(:), 'type', 'spearman');

rp = corr(xb, pb);
rp0 = corr(xb(:), pb(:));

ck = corr(xb, kb, 'type', 'spearman');
ck0 = corr(xb(:), kb(:), 'type', 'spearman');

rk = corr(xb, kb);
rk0 = corr(xb(:), kb(:));

fprintf('mode is %d\n', mode);
M = [[diag(cp)' cp0]; [diag(rp)' rp0]; [diag(ck)' ck0]; [diag(rk)' rk0]; ];
disp(num2cell(M))
%--------------------------------------------------------------------------
if nargin< 2
    close all;    
    nr = 1;
    nc = 2;
    fsiz = [0 0 .65 .25];
    subplots = 1:2;
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);    
end

fsy = def('fsy');
fs = def('fs');
cols = def('col');

% conditions = {[1 3], [2 4]};
% for k=1:2
%     subplot(nr, nc, subplots(k));
%     for j=1:2
%         i = conditions{k}(j);
%         plot_ci(xb(:, i), xe(:, i), cols(j, :));
%         ylim([.5 1])    
%         hold on;
%     end
% end

ii = {[1 3], [2 4]};
label = 'Learning rate coefficient';
lg_labels = {'Small', 'Large'};
lg_title = def('sto');
fig_titles = {sprintf('Small %s', lower(def('vol'))), sprintf('Large %s', lower(def('vol')))};
yl = [0.5 1];
for i=1:2
    subplot(nr, nc, subplots(i));
    for j=1:2
    hm(j) = plot_ci(xb(:, ii{i}(j)), xe(:, ii{i}(j)), cols(j, :)); hold on;
    ylim(yl)
    end
    ylabel(label, 'fontsize', fsy);
    xlabel('Trial', 'fontsize', fsy);
    title(fig_titles{i}, 'fontsize', fsy, 'fontweight', 'normal');
    set(gca,'YTick', .5:.1:1);

    if i==1
        lg = legend(hm, lg_labels, 'box', 'off', 'location', 'north', 'fontsize', fs);
        title(lg, lg_title,'fontweight','normal', 'fontsize', fsy);
    end
end


end

function dat = sim_model(data_set)
model_name = 'pf_Fd1w10_seed0';

if data_set == 1
    subjs = 1:223;
elseif data_set == 2
    subjs = 224:643;
elseif data_set == 0
    subjs = 1:643;
end

fitdir = fullfile('..', 'fit_pf', model_name);
% 
% fsig = fullfile(fitdir, sprintf('rev_pf_med_dynamics_data%d.mat', data_set));
fsig = fullfile(fitdir, sprintf('rev_pf4_dynamics_data%d.mat', data_set));


data = get_data(data_set);

tt = 1:50;
if 1 %~exist(fsig,'file')
    fnames = getfileordered(fitdir,'subj_%04d.mat',subjs);    
    N = length(fnames);    
    for n=1:N
        fname  = fnames{n};
        f = load(fname);
        parameters(n,:) = f.pf_observed.x;

        val = f.pf_observed.signal.dynamics_mean.val;
        bucket = data{n}.bucket;
        s = std(val-bucket);
        val = val + s.*randn(size(val));
        lr(:, :, n) = f.pf_observed.signal.dynamics_mean.lr;  
%         for k=1:10
%             sig(:, :, k) = f.pf_observed.signal.dynamics.val{k};
%         end
%         val = median(sig, 3);

        outcome = data{n}.bag;
        dat{n} = struct('bucket', val, 'bag', outcome);        
    end    
    
end
end

% % % % % % % 
function dat = sim_kalman(data_set)


[data] = get_data(data_set);

% load Kalman with 1 parameter
% fitdir = fullfile('..','analysis','fit_kf2');
% fname = fullfile(fitdir, sprintf('dataset%d_vol.mat', data_set));

fitdir = fullfile('..','fit_kpf','kalman4'); makedir(fitdir);
fname = fullfile(fitdir, sprintf('kalman_dataset%d.mat', data_set));
fsig = fullfile(fitdir, sprintf('rev_dynamics_data%d.mat', data_set));

f = load(fname);

lambda = f.kalman_parameter;

if 1 %~exist(fsig, 'file')
    for n=1:length(data)
        o = data{n}.bag;
        [val] = kalman(lambda(n,:), o);
        bucket = data{n}.bucket;
        s = std(val-bucket);
        val = val + s.*randn(size(val));        
        dat{n} = struct('bucket', val, 'bag', o);        
    end
end

end

function [val] = kalman(lambda, outcome)

ndim = size(outcome,2);
N = size(outcome,1);
lr = nan(N,ndim);

m = 60+zeros(1,ndim);
r = 1*ones(1,ndim);

for t=1:N
    val(t, :) = m;

    a = (r+lambda)./(r+1+lambda);
    m = m + a.*(outcome(t, :)-m);
    r = (1-a).*(r+lambda); 

    lr(t, :) = a;
end

c = corr(lr(1:end-1,:), lr(2:end, :));
c = diag(c);
end

