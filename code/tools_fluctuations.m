function [b, e] = tools_fluctuations(data)
T = 50;
w = 3;

if mod(w,2)==1
    t0 = ceil(w/2);
    l = floor(w/2);
    l_up = l;
else
    t0 = ceil(w/2) + 1;
    l = floor(w/2);
    l_up = 0;
end

for n=1:length(data)
    for t=1:T
        if t<t0
            tt = 1:(t+l);
        elseif t>(T-t0+1)
            tt = t-l:T;
        else
            tt = t-l:(t+l_up);
        end


        delta_all = []; update_all = [];     
        for j=1:4    
            dat = data{n};
            bucket = dat.bucket(tt,j);
            bag = dat.bag(tt,j);
            
            update = bucket(2:end)-bucket(1:end-1); % bucket(t) is an index of m(t-1) after the update
            delta = bag - bucket;
            delta = delta(1:end-1);
            
            delta_all = blkdiag(delta_all, delta);
            update_all = [update_all; update];
        end
        [bb] = glmfit(delta_all, update_all, 'normal','constant','off');

        a(t,:,n) = bb(1:end);
    end        
end
a(1:(t0-1),:,:) = [];

b = mean(a, 3);
e = serr(a, 3);

% b = median(a,3);
% for i=1:4
%     x = a(:,i,:);
%     siz = size(x);
%     x = reshape(x, siz([1 3]));
%     e(:, i) = se_median(x', 500);
% end
end

