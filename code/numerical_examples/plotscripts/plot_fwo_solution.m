prepare
load([plotdatapath 'fwo_solution.mat'])

f1 = setfig('b1');
hold off
scontour(xi/a, yi/a, svmi/10^6, 1000, 1000, 20, 'LineStyle', 'none')
box on
xlim([-3 3])
ylim([-3,0])
xlabel('$x/a$')
ylabel('$y/a$')
ch = colorbar;
title('von Misesova napetost')
title(ch, 'MPa', 'interpreter', 'latex');
daspect([1 1 1])

exportfig(f1, [imagepath 'fwo_solution.png'], '-png', '-r300')