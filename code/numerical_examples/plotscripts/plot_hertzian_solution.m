prepare
load([plotdatapath 'hertzian_solution.mat'])

f2 = setfig('b1');
plot(xt/b, sxxt/p0, 'x-');
plot(xt/b, asxxt/p0, '-');
xlabel('$x/b$')
ylabel('$y/b$')
legend('$\sigma_{xx}$ izra\v{c}unan', '$\sigma_{xx}$ analiti\v{c}en');

zoomx = [-1.13, -0.98];
zoomy = [-2e-3, 3e-3];
axpos = [.521 .16 .376 .399];
ax2 = zoomin(f2,...
       [zoomx(1) zoomy(1)-0.015 zoomx(2)-zoomx(1) zoomy(2)-zoomy(1)+0.03], ...
       axpos);

plot(ax2, xt/b, sxxt/p0, 'x-');
plot(ax2, xt/b, asxxt/p0, '-');
ylim(zoomy)


f3 = setfig('b3');
vmi = von_mises(sxxi, syyi, sxyi);

f = 5e3;
scatter((xi+f*ui)/b, (yi+f*vi)/b, 5, vmi/p0, 'filled')
colorbar
xlim([-3.2, 3.2])
daspect([1 1 1])
xlabel('$x/b$')
ylabel('$y/b$')
title('von Misesova napetost$/p_0$', 'interpreter', 'latex')


% plot(xt/b, syyt/p0, '.-');
% plot(xt/b, sxyt/p0, '.-');

% f1 = setfig('b2', 'Position', [2000, 500, 1000, 500]);
% hold off
% subplot(2, 2, 1)
% scontour(xi/b, yi/b, sxxi/p0, 100, 100, 30, 'LineStyle', 'none')
% xlabel('$x/b$')
% ylabel('$y/b$')
% title('$s_{xx}/p_0$', 'interpreter', 'latex')
% colorbar
% ax = gca;
% ax.Position(1) = 0.075;
% ax.Position(3) = 0.35;
% ax.Position(2) = 0.61;
% 
% subplot(2, 2, 2)
% scontour(xi/b, yi/b, sxyi/p0, 100, 100, 30, 'LineStyle', 'none')
% colorbar
% xlabel('$x/b$')
% ylabel('$y/b$')
% title('$s_{xy}/p_0$', 'interpreter', 'latex')
% ax = gca;
% ax.Position(3) = 0.35;
% ax.Position(1) = 0.56;
% ax.Position(2) = 0.61;
% 
% subplot(2, 2, 3)
% scontour(xi/b, yi/b, syyi/p0, 100, 100, 30, 'LineStyle', 'none')
% colorbar
% xlabel('$x/b$')
% ylabel('$y/b$')
% title('$s_{yy}/p_0$', 'interpreter', 'latex')
% ax = gca;
% ax.Position(1) = 0.075;
% ax.Position(3) = 0.35;
% 
% subplot(2, 2, 4)
% scontour(xi/b, yi/b, vmi/p0, 100, 100, 30, 'LineStyle', 'none')
% colorbar
% xlabel('$x/b$')
% ylabel('$y/b$')
% title('von Misesova napetost$/p_0$', 'interpreter', 'latex')
% ax = gca;
% ax.Position(3) = 0.35;
% ax.Position(1) = 0.56;

pause(15)
exportfig(f2, [imagepath 'hertzian_solution_sxx_top'], '-pdf')
pause(5)
exportfig(f3, [imagepath 'hertzian_solution_deformed_vm'], '-png', '-r300')