prepare
load([plotdatapath 'poisson_weird.mat'])

x = pos(1, :);
y = pos(2, :);
z = pos(3, :);

maxcol = max(sol);
mincol = min(sol);
colrat = (1 - (sol - mincol)/(maxcol-mincol)) / 3 * 2;
colors = hsv2rgb([colrat ones(size(colrat)) ones(size(colrat))]);


f1 = setfig('b1');
hold off
scatter3(x, y, z, 25, colors, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5)
daspect([1 1 1])
box off
axis off
az = -67.1000;
el = 22.8000;
view(az, el)
set(gcf,'Position', [0, 0, 1000 800])

exportfig(f1, [imagepath 'poisson_weird3'], '-png')