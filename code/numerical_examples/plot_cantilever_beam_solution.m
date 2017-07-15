prepare
load([plotdatapath 'cantilever_beam_solution.mat'])

f1 = setfig('b1');
subplot(2, 1, 1)
hold on
f = 1e5;
plot([0 L L 0 0], [-D/2 -D/2 D/2 D/2 -D/2], '--k')
scatter(x+f*u, y+f*v, 5, sxx/1000, 'filled')
% scontour(x+f*u, y+f*v, sxx/1000, 100, 100, 20);
c = colorbar;
title(c, 'kPa', 'interpreter', 'latex');
title('$\sigma_{xx}$')
daspect([1, 1, 1])
box on
grid on
ylim([-4.5 2.75])
xlim([-0.5 30.25])

subplot(2, 1, 2)
hold on
plot([0 L L 0 0], [-D/2 -D/2 D/2 D/2 -D/2], '--k')
scatter(x+f*u, y+f*v, 5, sxy/1000, 'filled')
% scontour(x+f*u, y+f*v, sxx/1000, 100, 100, 20);
c = colorbar;
title(c, 'kPa', 'interpreter', 'latex');
title('$\sigma_{xy}$')
daspect([1, 1, 1])
box on
grid on
ylim([-4.5 2.75])
xlim([-0.5 30.25])

exportfig(f1, [imagepath 'cantilever_beam_solution'], '-png', '-r300');