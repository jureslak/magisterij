prepare
datafile = [datapath 'fwo_table_wip.h5'];

name = '/rad0.01/cof0.3/calc0032';

pos = h5read(datafile, [name '/pos']);
stress = h5read(datafile, [name '/stress']);
N = h5readatt(datafile, name, 'N');

x = pos(1, :);
y = pos(2, :);

sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);
sv = von_mises(sxx, syy, sxy);

% displ  = h5read(datafile, [name '/disp']);
% u = displ(1, :);
% v = displ(2, :);
    
% max(max(displ))

a = h5readatt(datafile, name, 'a');
p0 = h5readatt(datafile, name, 'p0');

R = h5readatt(datafile, name, 'R');
COF = h5readatt(datafile, name, 'COF');

femnizname = [plotdatapath sprintf('fem_niz/surface_sxx_r%d_mu%g.txt', 1000*R, COF)];
femnizdata = dlmread(femnizname, '\t');

% M = spconvert(h5read(datafile, [name '/matrix'])');
% rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% max(max(displ2))

%%
close all

f1 = setfig('b1');
scatter(x/a, y/a, 5, sxx, 'filled')
daspect([1 1 1])

f2 = setfig('b2');
top = find(y == max(y));
[~, I] = sort(x(top));
top = top(I);

plot(x(top)/a, sxx(top)/10^6, 'o');
plot(femnizdata(:, 1)/a/10^3, femnizdata(:, 2), '.-');
xlim([-3, 3])



% ylim([-3, 3])