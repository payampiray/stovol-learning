function stats_table

clc;
over_write = 0;
experiment = 1;

i = 1;

st =  get_demo();
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@model_neutral, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = model_neutral_1st_block;
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
disp(T);

st = get_stats(@model_neutral_performance, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@kf_maladaptivity, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@hpl_parameters, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@other_bmc, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@other_parameters, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@hpl_modulation, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@hpl_fluctuations, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

% % % ---------------------------------------------------------------------


experiment = 2;
over_write = 0;

st = get_stats(@model_neutral, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@model_neutral_performance, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@kf_maladaptivity, experiment, over_write);
copy_table(st.table.data, 2);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

st = get_stats(@hpl_parameters, experiment, over_write);
copy_table(st.table.data, 3);
T = array2table(st.table.data, 'VariableNames', st.table.columns, 'RowNames', st.table.rows);
fprintf('\nSupplementary Table %d\n', i); i = i+1; disp(T);

end

function str = copy_table(x, n)
if length(n)==1, n = n*ones(size(x, 1), 1); end

for i=1:size(x, 1)
    y(i, :) = round(x(i, :)*10^n(i))/10^n(i);
end

str = num2clip(y);

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