function plot_raincloud(mb, eb, b, dens_norm)
if nargin<4, dens_norm = 1; end

fsy = def('fs');

colmap = def('col_unique');
config.dens_ratio = 1.5/dens_norm;
config.face_color = colmap;
config.patch_alpha = .4;
config.dot_displacement = .1;
config.dots_width = .1;
% config.bandwidth = .0001;
x_mid = 1;
x_bar = 0.6;
xl = [0 2];
barwidth = .15;

raincloud1xN(b(:, 1), x_mid, config)
set(gca, 'fontsize', fsy,'box', 'off', 'xlim', xl, 'xtick', [],'box', 'off');

hold on;
bar(x_bar, mb, barwidth,'FaceColor', colmap(1,:),'EdgeColor','k','linewidth',1, 'FaceAlpha',1); hold on;
plot(x_bar*ones(1,2), mb+[-eb eb], 'color','k','linewidth',2); hold on;
end
