prepare
datafile = [datapath 'hertzian_refined_convergence_jarjar.h5'];
info = h5info(datafile);

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');

grp = info.Groups(7).Groups(19);
name = grp.Name;

N = h5readatt(datafile, name, 'N');
pos = h5read(datafile, [name '/pos']);
x = pos(1, :);
y = pos(2, :);

[asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);

stress = h5read(datafile, [name '/stress']);
sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);

displ = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);

testidx = find(x.^2 + y.^2 < (1000*b)^2);% & y.^2 > (b/2)^2);

errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
erryy = max(max(abs(syy(testidx) - asyy(testidx))));
errxy = max(max(abs(sxy(testidx) - asxy(testidx))));

err = max([errxx, erryy, errxy])/p0

lim = 3*b;
interesting = -lim < x & x < lim & -lim < y;

xi = x(interesting);
yi = y(interesting);
sxxi = sxx(interesting);
syyi = syy(interesting);
sxyi = sxy(interesting);
ui = u(interesting);
vi = v(interesting);

top = find(y == max(y) & -lim < x & x < lim);
[~, I] = sort(x(top));
top = top(I);

xt = x(top);
yt = y(top);
sxxt = sxx(top);
syyt = syy(top);
sxyt = sxy(top);
[asxxt, asyyt, asxyt] = hertzian_analytical(x(top), y(top), b, p0); 

save([plotdatapath 'hertzian_solution.mat'], 'b', 'p0',...
    'xi', 'yi', 'sxxi', 'sxyi', 'syyi', 'ui', 'vi',...
    'xt', 'yt', 'sxxt', 'sxyt', 'syyt', 'asxxt');