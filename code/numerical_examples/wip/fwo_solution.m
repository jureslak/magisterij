prepare
datafile = [datapath 'fwo_solution_wip.h5'];
info = h5info(datafile);

a = h5readatt(datafile, '/', 'a');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');

name = info.Groups(1).Groups(1).Name;

pos = h5read(datafile, [name '/pos']);
x = pos(1, :);
y = pos(2, :);

stress = h5read(datafile, [name '/stress']);
sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);
svm = von_mises(sxx, syy, sxy);

displ  = h5read(datafile, [name '/displacement']);
u = displ(1, :);
v = displ(2, :);


lim = 3*a;
interesting = -lim < x & x < lim & -lim < y;
xi = x(interesting);
yi = y(interesting);
sxxi = sxx(interesting);
syyi = syy(interesting);
sxyi = sxy(interesting);
svmi = svm(interesting);
ui = u(interesting);
vi = v(interesting);

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


f2 = setfig('b2');
f = 1e2;
scatter((xi+f*ui)/a, (yi+f*vi)/a, 5, svmi, 'filled')
% xlim([-3 3])
% ylim([-3,0])
xlabel('$x/a$')
ylabel('$y/a$')
ch = colorbar;
title('von Misesova napetost')
title(ch, 'MPa', 'interpreter', 'latex');
daspect([1 1 1])

% exportfig(f1, '../../../images/hertzian_convergence', '-pdf')
% exportfig(f2, '../../../images/hertzian_time', '-pdf')