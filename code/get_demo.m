function st = get_demo()

[data] = get_data(1);
for n=1:length(data)
    demo.age(n,1) = data{n}.age;
    demo.gender(n,1) = data{n}.gender;
end

n1 = length(demo.age);
mx1 = nanmean(demo.age);
sx1 = nanstd(demo.age);
nx1 = nansum(demo.gender==1);
women1 = nansum(demo.gender=='F');
men1 = nansum(demo.gender=='M');
unreport1 = length(data) - women1 - men1;


[data] = get_data(2);
for n=1:length(data)
    demo.age(n,1) = data{n}.age;
    demo.gender(n,1) = data{n}.gender;
end

n2 = length(data);
mx2 = nanmean(demo.age);
sx2 = nanstd(demo.age);
women2 = nansum(demo.gender=='F');
men2 = nansum(demo.gender=='M');
unreport2 = length(data) - women2 - men2;

st.table.data = [[n1; mx1; sx1; women1; men1; unreport1] [n2; mx2; sx2; women2; men2; unreport2]];

st.table.columns = {'Experiment 1','Experiment 2'};
st.table.rows = {
'Sample size'
'Mean age'
'SD age'
'Gender: female'
'Gender: male'
'Gender: none or unreported'};

end