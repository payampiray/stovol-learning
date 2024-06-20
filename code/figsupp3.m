function figsupp3

h = other_signals(1);

% --------
fn = def('fn');
fsA = def('fsA')+2;
xsA = -.25;def('xsA');
ysA = 1.17;

abc = 'abcd';
h = h(1:4:end);

for i= 1:length(h)
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i),'fontweight','bold');
end
end