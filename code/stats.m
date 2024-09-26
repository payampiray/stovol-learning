function stats
clc;
over_write = 0;
experiment = 1;

fprintf('%s\n\n', repmat('=',1,40));
fprintf('Experiment 1\n');
fprintf('%s\n\n', repmat('=',1,40));

st = get_stats(@model_neutral_performance, experiment, over_write);
write_it('performance', st, 1:3)

st = get_stats(@model_neutral_performance_blk, experiment, over_write);
write_it('performance relative to var', st, 1:4)

st_glm = get_stats(@model_neutral, experiment, over_write);
write_it('model neutral', st_glm, 2:4);

st = get_stats(@model_neutral_1st_block, experiment, over_write);
write_it('model_neutral_1st_block', st, 1:10)

st_bmc = get_stats(@kf_bmc, experiment, over_write);
write_it('BMC between KF and RW', st_bmc, 1:2);

st_kalman = get_stats(@kf_maladaptivity, experiment, over_write);
write_it('Kalman: lambda_s', st_kalman.lambda_s.all, 1);
write_it('Kalman: lambda_v', st_kalman.lambda_v.all, 1);
write_it('The maladaptive group based on lambda_s', st_kalman.lambda_s.neg, 2)
write_it('The maladaptive group based on lambda_v', st_kalman.lambda_v.neg, 1)


st_bmc = get_stats(@hpl_bmc, experiment, over_write);
write_it('BMC between PF and KF', st_bmc, 1:2);

% write_it('model agnostic BF01 for the interaction', st_glm, 4);
% write_it('Kalman: lambda_i', st_kalman.lambda_i.all, 1);

st_fluctions = get_stats(@hpl_fluctuations, experiment, over_write);
write_it('correlation fluctuations', st_fluctions.all, 1);


st_signal = get_stats(@hpl_signal, experiment, over_write);
write_it('PF: lr estimate', st_signal.lr, 1:2);
write_it('PF: sto estimate', st_signal.sto, 1:2);
write_it('PF: vol estimate', st_signal.vol, 1:2);

st_dynamics = get_stats(@hpl_dynamics, experiment, over_write);
write_it('PF: vol dynamics', st_dynamics.vol, 3);
write_it('PF: sto dynamics', st_dynamics.sto, 3);
write_it('PF: difr dynamics', st_dynamics.difr, 3);

st_clustering = get_stats(@hpl_clustering, experiment, over_write);
write_it('Model clustering analysis', st_clustering.a_model, 1:2);
write_it('Data clustering analysis', st_clustering.a_data, 1:2);

st_assignment = get_stats(@hpl_assignment, experiment, over_write);
write_it('Model |LR| against |AC|', st_assignment.model, 1:2);
write_it('Model |DeltaV/DeltaS| against |AC|', st_assignment.b2v, 1:2);
write_it('Data |LR| against |AC|', st_assignment.data, 1:2);

st_rt = get_stats(@model_neutral_response_time, experiment, over_write);
write_it('RT against |AC|', st_rt.data, 1:2);

st_modulation = get_stats(@hpl_modulation, experiment, 1);
write_it('Model modulation analysis', st_modulation.model, 5:6);
write_it('Data modulation analysis', st_modulation.data, 5:6);

% -------------------------------------------------------------------------
fprintf('%s\n\n', repmat('=',1,40));
fprintf('Experiment 2\n');
fprintf('%s\n\n', repmat('=',1,40));

experiment = 2;


st_performance = get_stats(@model_neutral_performance, experiment, over_write);
write_it('performance', st_performance, 1:3)

st_glm = get_stats(@model_neutral, experiment, over_write);
write_it('model neutral', st_glm, 2:4);

st_kalman = get_stats(@kf_maladaptivity, experiment, over_write);
write_it('Kalman: lambda_s', st_kalman.lambda_s.all, 1)
write_it('Kalman: lambda_v', st_kalman.lambda_v.all, 1)
write_it('The maladaptive group based on lambda_s', st_kalman.lambda_s.neg, 2)
write_it('The maladaptive group based on lambda_v', st_kalman.lambda_v.neg, 1)

st_signal = get_stats(@hpl_signal, experiment, over_write);
write_it('PF: sto estimate', st_signal.sto, 1:2);
write_it('PF: vol estimate', st_signal.vol, 1:2);

st_clustering = get_stats(@hpl_clustering, experiment, over_write);
write_it('Model clustering analysis', st_clustering.a_model, 1:2);
write_it('Data clustering analysis', st_clustering.a_data, 1:2);

st_assignment = get_stats(@hpl_assignment, experiment, 1);
write_it('Model |LR| against |AC|', st_assignment.model, 1:2);
write_it('Model |DeltaV/DeltaS| against |AC|', st_assignment.b2v, 1:2);
write_it('Data |LR| against |AC|', st_assignment.data, 1:2);

st_rt = get_stats(@model_neutral_response_time, experiment, 1);
write_it('RT against |AC|', st_rt.data, 1:2);

st_modulation = get_stats(@hpl_modulation, experiment, 1);
write_it('Model modulation analysis', st_modulation.model, 5:6);
write_it('Data modulation analysis', st_modulation.data, 5:6);


end

function write_it(st_name, stats, idx)
fnames = fieldnames(stats);
fnames(strcmp(fnames, 'table')) = [];

fprintf('%s\n', st_name);
for i=1:length(fnames)
    x= stats.(fnames{i});
    if size(x,1)==1
        field = fnames{i}; 
        st.(field) = x(idx);
    else
        for j=1:size(x, 1)
            field = sprintf('%s_%d', fnames{i}, j);
            st.(field) = x(j, idx);
        end
    end
    
end
display(st);
fprintf('%s\n\n', repmat('-',1,40));
end

function stats = get_stats(analysis_mfile, experiment, over_write)
if nargin<3
    over_write = 0;
end

mfile = func2str(analysis_mfile);
fname = fullfile(sprintf('experiment%d', experiment), sprintf('%s.mat', mfile));
if ~exist(fname, 'file') || over_write
    stats = analysis_mfile(experiment); close all;
    save(fname, 'stats');    
end
fst = load(fname);
stats = fst.stats;
end
