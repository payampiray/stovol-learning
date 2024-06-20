function h = raincloud1xN(y, y_base, config)
K = size(y,2);
if nargin<3
    config = [];
end
if nargin<4
    is_right = ones(1,K);
end


if isempty(config)
    config = struct('K', K);
end

p = inputParser;
p.addParameter('K', K);
p.addParameter('is_right', ones(1,K));

p.addParameter('support', 'unbounded');
p.addParameter('bandwidth', nan);
p.addParameter('n_points', 50);
p.addParameter('dens_ratio', 0.25);
p.addParameter('n_extra_y', 0.04);
p.addParameter('patch_alpha', 0.6);
p.addParameter('face_color', 0.3*ones(K,3));
p.addParameter('add_patch', 1);


p.addParameter('dots_width', 0.2);
p.addParameter('dot_displacement', 0.15);
p.addParameter('dots_size', 5);
p.addParameter('dots_color', 0.4*ones(1,3));

p.addParameter('xlim_margin', 0.5);


p.parse(config);
config    = p.Results;

for k=1:K
    ax = y_base;    
    h(k) = draw(ax, y(:,k), is_right(k), config.face_color(k,:), config);
    hold on;
end
% % 
% % set(gca, 'xtick', 1:ax);
% % xlim([1 ax] + [-1 1]*config.xlim_margin);
% % set(gca,'xticklabel',labels);

end

function h = draw(ax, y, is_right, face_color, config)
box_mid = ax;


% patch
add_patch = config.add_patch;
support = config.support;
bandwidth = config.bandwidth;
dens_ratio = config.dens_ratio;
n_extra_y = config.n_extra_y;
n_points = config.n_points;
patch_alpha = config.patch_alpha;
 
% % bar
% bar_mid = .05;
% bar_width = .3; % width
% bar_linewidth = 1; % line width
% rect_alpha = alf;

% dots
dots_width = config.dots_width;
dot_displacement = config.dot_displacement;
dots_size = config.dots_size;
dots_color = config.dots_color;
%------------------------------------
is_right(is_right==0) = -1;

N = length(y);
maxy=max(y);
miny=min(y);

extra_y = n_extra_y*(maxy-miny)/(n_points);
binranges = linspace(miny-extra_y,maxy+extra_y,n_points);

if isnan(bandwidth)
    [dens,dens_pos, bw] = ksdensity(y,binranges,'function','pdf','kernel','normal', 'support', support); %#ok<ASGLU> 
else
    [dens,dens_pos] = ksdensity(y,binranges,'function','pdf','kernel','normal','Bandwidth', bandwidth, 'support', support);
end
area=sum(dens)*2*(binranges(2)-binranges(1));

%Normalized area is a third of the max available area
max_area= maxy-miny;
dens=dens*max_area/(3*area)*dens_ratio;
% plot(-dens,1:npoints);

% if left
if is_right == -1
    xpatch=[box_mid-dens(1:end-1) ; ...
            box_mid+zeros(1,length(dens)-1) ; ...
            box_mid+zeros(1,length(dens)-1);  ...
            box_mid-dens(2:end)];
end

% if right    
if is_right == 1
    xpatch=[box_mid+zeros(1,length(dens)-1) ; ...
            box_mid+dens(1:end-1) ; ...
            box_mid+dens(2:end) ; ...
            box_mid+zeros(1,length(dens)-1)];    
end

ypatch=[dens_pos(1:end-1) ; ...
    dens_pos(1:end-1) ; ...
    dens_pos(2:end) ; ...
    dens_pos(2:end)];

minxp = min(min(xpatch,[],2));
maxxp = max(max(xpatch,[],2));
% rxp = ceil((maxxp - minxp)*10)/10;

if add_patch
hpatch = patch(xpatch,ypatch,[1 1 1],'EdgeAlpha',1,'FaceColor',face_color,'EdgeColor','none','FaceAlpha',patch_alpha);    
hpatchl = line(box_mid +is_right*[0 dens 0 0], [dens_pos(1) dens_pos dens_pos(end) dens_pos(1)],'color',zeros(1,3));
% xl1 = rxp*[-1 1]+box_mid;
hold on;
end

% bar_mid = box_mid + -isright*(bar_mid + bar_width);
% hrect = rectangle('position',[bar_mid-bar_width/2 0 bar_width mdy],'FaceColor',[face_color rect_alpha],'linewidth',bar_linewidth);
% hrectl = plot([bar_mid bar_mid],mdy+[-edy edy],'-','color','k','linewidth',2);

dots_mid = box_mid -is_right*dot_displacement;
xp = (rand(1,N)*dots_width - dots_width/2) +dots_mid;
hscat = scatter(xp,y,dots_size,repmat(dots_color,N,1),'filled');

h = gca;

xtick = get(gca,'xtick');
end