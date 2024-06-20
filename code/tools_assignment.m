function [b, z, bv2s] = tools_assignment(action, outcome, vol, sto)

dlr_all = [];
acov = [];
for j=1:4
    bucket = action(:,j);
    bag = outcome(:,j);
    
    delta = bag - bucket;
    delta = delta(1:end-1);

    acov = cat(1, acov, delta.*[0; delta(1:end-1)]);
    
    update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
    delta = bag - bucket;
    delta = delta(1:end-1);
    lr = update./delta;
    lr(lr<0) = 0;
    lr(lr>1) = 1;
 
    dlr = lr(2:end) - lr(1:end-1);
    dlr_all = cat(1, dlr_all, [nan; dlr]);
end

acov_abs = abs(acov);

a = glmfit([acov_abs acov] , abs(dlr_all) );
b = [a(2:end); a(1)];

% xt = prctile(acorr_abs, 10:10:90);
xt = prctile(acov_abs, 20:20:80);
y = abs(dlr_all);

xt = [0 xt max(acov_abs)];
for i=2:(length(xt))
    t = (acov_abs>=xt(i-1)) & (acov_abs<xt(i));
    z(i-1) = nanmean(y(t));
end

bv2s = [];
if nargin>2
    dvol_all = [];
    dsto_all = [];
    vol(end, :) = [];
    sto(end, :) = [];
    for j=1:4
        dvol = vol(2:end, j) - vol(1:end-1, j);
        dvol = [0;dvol];
        dvol_all = [dvol_all; dvol];
        dsto = sto(2:end, j) - sto(1:end-1, j);
        dsto = [dsto; 0];
        dsto_all = [dsto_all; dsto];

    end

    a = glmfit([acov_abs acov] , abs(log(dvol_all+eps)-log(dsto_all+eps)) );
    bv2s = [a(2:end); a(1)];
end

end
