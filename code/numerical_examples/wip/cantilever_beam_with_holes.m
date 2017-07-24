prepare
datafile = [datapath 'cantilever_beam_with_holes_wip.h5'];

info = h5info(datafile);
name = info.Groups(1).Groups(1).Name;

pos = h5read(datafile, [name '/pos']);
types = h5read(datafile, [name '/types']);
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

displ  = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);
max(sqrt(u.^2+v.^2))
 
    
M = spconvert(h5read(datafile, [name '/matrix'])');
rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% um = displ2(1:N);
% vm = displ2(N+1:end);
% max(sqrt(um.^2+vm.^2))

size(M)
%%
close all

sf1 = setfig('b1');
f = 1e5;
% plot([0 L L 0 0], [-D/2 -D/2 D/2 D/2 -D/2], '--k')
% c = zeros(size(x));
% c(types < 0) = 1;
% c(types == 1) = -1;

sm = sum(abs(M),2);
% types(1) = 3;
% scatter(x, y, 5, types, 'filled')
scatter(x+f*u, y+f*v, 5, sv/1000, 'filled')
colormap parula
% brighten(-0.7)
% scontour(x+f*u, y+f*v, sxx/1000, 100, 100, 50);
c = colorbar;
title(c, 'kPa', 'interpreter', 'latex');
title('$\sigma_{xx}$')
daspect([1, 1, 1])
ylim([-4.5 2.75])
xlim([-0.5 30.25])

