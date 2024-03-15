function [p, H, st] = signrank2(x)
n = size(x,2);

p = nan(1,n);
H = nan(1,n);
for i=1:n
    [p(i),H(i),stats(i)] = signrank(x(:,i)); %#ok<AGROW>
end

st = stats(1);
ffn = fieldnames(st);
for i=2:n
    for j=1:length(ffn)
        st.(ffn{j}) = [st.(ffn{j}) stats(i).(ffn{j})];
    end
end

end