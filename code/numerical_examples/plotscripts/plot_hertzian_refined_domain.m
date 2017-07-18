prepare
load([plotdatapath 'hertzian_refined_domain.mat']);

f1 = setfig('b1', 'Position', [2000, 1000, 600, 300], 'Visible', 'on');
gray = [1 1 1]*0.4;

scatter(xi/b, yi/b, 1, gray, 'filled')
caxis([-1e-3, 1e-3])
xlabel('$x/b$')
ylabel('$y/b$')

% daspect([1 1 1])
axpos = [.52 .18 .376 .399]; % if colorbar present change first to .47
ax = zoomin(f1, [-1.25 -0.25 0.5 0.25], axpos);
scatter(xi/b, yi/b, 1, gray, 'filled')

exportfig(f1, [imagepath 'hertzian_refined_domain'], '-png', '-r300')