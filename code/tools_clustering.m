function [mx, ms, mv, labels] = tools_clustering(action, outcome, sto, vol)

for j=1:4
    
    bucket = action(:,j);
    bag = outcome(:,j);
    
    update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
    delta = bag - bucket;
    delta = delta(1:end-1);

    a = update./delta;
    a(a<0) = 0;
    a(a>1) = 1;    
    lr(:, j) = a;
    PE(:, j) = [0; (delta(1:end-1))];
end

dlr = lr(2:end, :) - lr(1:end-1, :);
PE(end, :) = [];

sample_autocorr = PE(2:end, :).*PE(1:end-1, :);
dlr(end, :) = [];

t_neg = sample_autocorr<0;
t_pos = sample_autocorr>0;

x = dlr;
mx = nan(2, 4);
for i=1:4
    mx(1, i) = mean(x(t_neg(:, i), i));
    mx(2, i) = mean(x(t_pos(:, i), i));
end
mx = mean(mx, 2)';

labels = {'Negative', 'Positive'};

ms = [];
mv = [];
if nargin>2
    sto(end, :) = [];
    vol(end, :) = [];
    dsto = sto(2:end, :) - sto(1:end-1, :);
    dvol = vol(2:end, :) - vol(1:end-1, :);
    dsto(end, :) = [];
    dvol(end, :) = [];
        
    xx = {dsto, dvol};
    for j=1:2
        x = xx{j};
        for i= 1:4
            m_signal{j}(1, i) = mean(x(t_neg(:, i), i));
            m_signal{j}(2, i) = mean(x(t_pos(:, i), i));
        end
    end
ms = mean(m_signal{1}, 2)';
mv = mean(m_signal{2}, 2)';
end

end