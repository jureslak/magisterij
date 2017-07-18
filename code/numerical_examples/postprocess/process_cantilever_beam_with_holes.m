prepare
datafile = [datapath 'cantilever_beam_with_holes.h5'];

info = h5info(datafile);
name = info.Groups(1).Groups(1).Name;

pos = h5read(datafile, [name '/pos']);
stress = h5read(datafile, [name '/stress']);
N = h5readatt(datafile, name, 'N');

P = h5readatt(datafile, '/', 'P');
D = h5readatt(datafile, '/', 'D');
L = h5readatt(datafile, '/', 'L');
E = h5readatt(datafile, '/', 'E');
nu = h5readatt(datafile, '/', 'v');

lam = h5readatt(datafile, '/', 'lambda');
mu = h5readatt(datafile, '/', 'mu');

x = pos(1, :);
y = pos(2, :);

sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);
sv = von_mises(sxx, syy, sxy);
sh = tau1(sxx, syy, sxy);

displ  = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);
max(sqrt(u.^2+v.^2))

save([plotdatapath 'cantilever_beam_with_holes.mat'], 'x', 'u', 'y', 'v',...
     'sv', 'L', 'D');

