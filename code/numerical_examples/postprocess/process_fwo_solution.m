prepare
datafile = [datapath 'fwo_solution.h5'];
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

lim = 4*a;
interesting = -lim < x & x < lim & -lim < y;
xi = x(interesting);
yi = y(interesting);
sxxi = sxx(interesting);
syyi = syy(interesting);
sxyi = sxy(interesting);
svmi = svm(interesting);
ui = u(interesting);
vi = v(interesting);

save([plotdatapath 'fwo_solution.mat'], 'xi', 'yi', 'svmi', 'a');