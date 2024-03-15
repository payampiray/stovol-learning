function st = other_parameters(experiment)
if nargin<1, experiment = 1; end

tbl = other_rw(experiment);

models = { @other_dbd1, @other_dbd4, @other_weber, @other_hgf};
for i=1:length(models)
    tbl_model = models{i}(experiment);
    tbl.data = [tbl.data; tbl_model.data];
    tbl.rows = [tbl.rows tbl_model.rows];
end

st.table = tbl;
end