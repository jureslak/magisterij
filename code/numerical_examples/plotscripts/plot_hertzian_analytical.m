prepare
% no input data

b = 1.3073e-04;
p0 = 2.6443e+06;
x = linspace(-3*b, 3*b, 100);
y = linspace(-3*b, 0, 100);
[X, Y] = meshgrid(x, y);
[sxx, syy, sxy] = hertzian_analytical(X, Y, b, p0);

close all
f1 = setfig('b1', 'Position', [2000, 500, 1000, 500]);
hold off
subplot(2, 2, 1)
contourf(X/b, Y/b, sxx/p0, 30, 'LineStyle', 'none')
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{xx}/p_0$', 'interpreter', 'latex')
colorbar
ax = gca;
ax.Position(1) = 0.075;
ax.Position(3) = 0.35;
ax.Position(2) = 0.61;

subplot(2, 2, 2)
contourf(X/b, Y/b, sxy/p0, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{xy}/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(3) = 0.35;
ax.Position(1) = 0.56;
ax.Position(2) = 0.61;

subplot(2, 2, 3)
contourf(X/b, Y/b, syy/p0, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('$s_{yy}/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(1) = 0.075;
ax.Position(3) = 0.35;

vm = von_mises(sxx, syy, sxy);
subplot(2, 2, 4)
contourf(X/b, Y/b, vm/p0, 30, 'LineStyle', 'none')
colorbar
xlabel('$x/b$')
ylabel('$y/b$')
title('von Mises stress$/p_0$', 'interpreter', 'latex')
ax = gca;
ax.Position(3) = 0.35;
ax.Position(1) = 0.56;

exportfig(f1, [imagepath 'hertzian_analytical'], '-opengl', '-png', '-r300')