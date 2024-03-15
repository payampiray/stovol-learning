function figsupp2

h = other_signals(1);

% --------
fn = def('fn');
fsA = def('fsA')+2;
xsA = -.2;def('xsA');
ysA = 1.1;

abc = 'abcdefghijkl';

for i= 1:length(h)
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i),'fontweight','bold');
end

end