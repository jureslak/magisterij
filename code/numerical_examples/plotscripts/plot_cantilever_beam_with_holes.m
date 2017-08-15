prepare
load([plotdatapath 'cantilever_beam_with_holes.mat']);
 
f1 = setfig('b1');
f = 1e5;
plot([0 L L 0 0], [-D/2 -D/2 D/2 D/2 -D/2], '--k')
scatter(x+f*u, y+f*v, 1, sv/1000, 'filled')
title('von Misesova napetost')
% scontour(x+f*u, y+f*v, sxx/1000, 100, 100, 50);
c = colorbar;
title(c, 'kPa', 'interpreter', 'latex');
daspect([1, 1, 1])
ylim([-4.5 2.75])
xlim([-0.5 30.25])

maxmove = max(sqrt(u.^2+v.^2))

exportfig(f1, [imagepath 'cantilever_beam_with_holes'], '-png', '-r300', '-opengl')