close all
clear
datafile = '../data/poisson_neumann.h5';

name = '/calc0101';

pos = h5read(datafile, [name '/pos']);
sol = h5read(datafile, [name '/sol']);
N = h5readatt(datafile, name, 'N');

x = pos(1, :);
y = pos(2, :);

setfig 'b1';
scatter3(x, y, sol, 25, sol, 'filled')
view([1, 1, 1])