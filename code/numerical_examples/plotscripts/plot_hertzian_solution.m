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
zoomxc = zoomx;
zoomyc = [-2e-3-0.015, 3e-3+0.015];
% small square
lc = [1 1 1]*0.3;
plot([zoomxc(1) zoomxc(1) zoomxc(2) zoomxc(2) zoomxc(1)],...
     [zoomyc(2) zoomyc(1) zoomyc(1) zoomyc(2) zoomyc(2)],...
     '-k', 'LineWidth', 0.5, 'Color', lc)

axpos = [.521 .16 .376 .399];
% lines
mrg = 0.002;
[xf, yf] = ds2nfu(zoomxc(1), zoomyc(1));
annotation('line', [xf axpos(1)], [yf axpos(2)+mrg], 'Linewidth', 0.5, 'Color', lc)
[xf, yf] = ds2nfu(zoomxc(2), zoomyc(2));
annotation('line', [xf axpos(1)+axpos(3)-mrg*0.5], [yf axpos(2)+axpos(4)-mrg], 'Linewidth', 0.5, 'Color', lc)
[xf, yf] = ds2nfu(zoomxc(1), zoomyc(2));
annotation('line', [xf axpos(1)], [yf axpos(2)+axpos(4)], 'Linewidth', 0.5, 'Color', lc)
[xf, yf] = ds2nfu(zoomxc(2), zoomyc(1));
annotation('line', [xf axpos(1)+axpos(3)-mrg], [yf axpos(2)+mrg], 'Linewidth', 0.5, 'Color', lc)
% annotation('rectangle', axpos, 'FaceColor', 'red', 'LineWidth', 0.2)

% 2nd plot
p = uipanel(f2, 'Position', axpos, 'BackgroundColor', 'white',...
            'HighLightColor', lc, 'BorderWidth', 1, 'BorderType', 'line');
ax2 = axes('Parent', p);
axis tight
box on % put box around new pair of axes
grid on
hold on
plot(ax2, xt/b, sxxt/p0, 'x-');
plot(ax2, xt/b, asxxt/p0, '-');
xlim(zoomx)
ylim(zoomy)

% plot(xt/b, syyt/p0, '.-');
% plot(xt/b, sxyt/p0, '.-');


f1 = setfig('b1', 'Position', [2000, 500, 1000, 500]);
hold off
subplot(2, 2, 1)
scontour(xi/b, yi/b, sxxi/p0, 100, 100, 30, 'LineStyle', 'none')
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{xx}/p_0$', 'interpreter', 'latex')
colorbar
ax = gca;
ax.Position(1) = 0.075;
ax.Position(3) = 0.35;
ax.Position(2) = 0.61;

subplot(2, 2, 2)
scontour(xi/b, yi/b, sxyi/p0, 100, 100, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{xy}/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(3) = 0.35;
ax.Position(1) = 0.56;
ax.Position(2) = 0.61;

subplot(2, 2, 3)
scontour(xi/b, yi/b, syyi/p0, 100, 100, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{yy}/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(1) = 0.075;
ax.Position(3) = 0.35;

subplot(2, 2, 4)
vmi = von_mises(sxxi, syyi, sxyi);
scontour(xi/b, yi/b, vmi/p0, 100, 100, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('von Misesova napetost$/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(3) = 0.35;
ax.Position(1) = 0.56;

% exportfig(f1, [imagepath 'hertzian_solution'], '-opengl', '-png', '-r100, 100, 300')

exportfig(f2, [imagepath 'hertzian_solution_sxx_top'], '-pdf')