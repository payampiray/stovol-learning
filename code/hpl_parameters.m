function st = hpl_parameters(experiment)
if nargin<1, experiment = 1; end

f = load(fullfile('..',sprintf('experiment%d', experiment),'model_hpl.mat'));
parameters = f.parameters.Variables;
x = prctile(parameters, [25 50 75]);

st.table.data = x;
st.table.columns = f.parameters.Properties.VariableNames;
st.table.rows = {'25% quantile', 'Median','75% quantile'};

end