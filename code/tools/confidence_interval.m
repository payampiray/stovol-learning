function ci = confidence_interval(x,dim)
if nargin<2, dim=1; end
[~, ~, ci] = ttest(x, 0, 'dim', dim);
ci = ci(2, :) - mean(x, dim);

end