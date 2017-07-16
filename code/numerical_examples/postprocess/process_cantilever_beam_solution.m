prepare
datafile = [datapath 'cantilever_convergence_jarjar.h5'];

name = '/mon9/calc0096';

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

displ  = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);
maxmove = max(sqrt(u.^2+v.^2)) % analitično izračunano v mathematici

[asxx, asyy, asxy, au, av] = cantilever_beam_analytical(x, y, P, L, D, E, nu);

errxx = max(max(abs(sxx - asxx)));
erryy = max(max(abs(syy - asyy)));
errxy = max(max(abs(sxy - asxy)));
Ms = max([max(abs(sxx)) max(abs(syy)) max(abs(sxy))]);

errstress =  max([errxx, erryy, errxy]) / Ms

Mu = max(sqrt(au.^2 + av.^2));
erru = max(max(abs(u - au)))/Mu;
errv = max(max(abs(v - av)))/Mu;
erruv = max([erru, errv])

save([plotdatapath 'cantilever_beam_solution.mat'], 'x', 'y', 'sxx', 'sxy', 'syy',...
    'u', 'v', 'N', 'maxmove', 'errstress', 'erruv', 'L', 'D');