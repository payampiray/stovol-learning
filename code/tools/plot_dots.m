function plot_dots(h, x, add_within_lines, col, dot_dist, marker_size)
% h: axes handle
% x: a cell with NxK matrix elements
% for a 2x2 plot, x{1} and x{2} should be Nx2

if nargin<3
    add_within_lines = false;
end
if nargin<4
    col = def('col');
end
if nargin<5
    dot_dist = .12;
end
if nargin<6
    marker_size = 7;
end

xt = get(h,'xtick');
xd = dot_dist*[-1 1];
dot_width = 0.03;


for j=1:length(xt)
    N = size(x{j},1);
    x_dots = xt(j)+xd+dot_width*randn(N,1);
    y_dots = x{j};
    if add_within_lines
        for n=1:N
            plot(h,x_dots(n,:),y_dots(n,:),'color',0.5*[1 1 1],'markersize',marker_size);
        end
    else
        for k=1:size(y_dots, 2)
            plot(h,x_dots(:,k),y_dots(:,k),'.','color',col(k,:),'markersize',marker_size);
        end
    end
end