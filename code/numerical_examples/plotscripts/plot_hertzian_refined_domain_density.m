prepare
load([plotdatapath 'hertzian_refined_domain_density.mat']);

f1 = setfig('b1', 'Position', [2000, 1000, 600, 300]);
hold off
xden = linspace(min(xi), max(xi), 500);
yden = linspace(min(yi), max(yi), 500);
[X,Y] = meshgrid(xden, yden);
level = log(avgd)/log(2);
Z = griddata(xi, yi, level, X, Y, 'linear');
contourf(X/b, Y/b, Z, 10, 'Linestyle', 'none');
caxis([min(level), max(level)]);
xlabel('$x/b$')
ylabel('$y/b$')
title('$\log_2(r_c/r_{c}^{\max})$')
colorbar
colormap(flipud(colormap))

exportfig(f1, [imagepath 'hertzian_refined_domain_density'], '-png', '-r300');