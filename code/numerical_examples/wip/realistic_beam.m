prepare
datafile = [datapath 'realistic_beam.h5'];

name = '/mon9/calc0100';

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

left = h5read(datafile, [name '/left'])+1;
right = h5read(datafile, [name '/right'])+1;
top = h5read(datafile, [name '/top'])+1;
bottom = h5read(datafile, [name '/bottom'])+1;
middle = h5read(datafile, [name '/middle'])+1;

x = pos(1, :);
y = pos(2, :);

sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);
sv = von_mises(sxx, syy, sxy);

displ  = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);
    
max(max(abs(displ)))

% M = spconvert(h5read(datafile, [name '/matrix'])');
% rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% msu = displ2(1:N);
% msv = displ2(N+1:end);

% errum = norm(msu'-u)
% errvm = norm(msv'-v)

% M = spconvert(h5read(datafile, [name '/matrix'])');
% rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% max(max(displ2))

%%
close all

f1 = setfig('b1');
c = zeros(size(x));
c(left) = 1;
c(right) = -1;
c(middle) = 3;
c(top) = 2;
c(bottom) = -2;

scatter(x, y, 5, c, 'filled')
colormap 'jet'
colorbar
daspect([1, 1, 1])
title('sxx')


f2 = setfig('b2');
scatter(x, y, 5, syy, 'filled')
colorbar
daspect([1, 1, 1])
title('syy')


f3 = setfig('b3');
scatter(x, y, 5, sxy, 'filled')
colorbar
daspect([1, 1, 1])
title('sxy')



f4 = setfig('b4');
scatter(x, y, 5, sxx, 'filled')
colorbar
daspect([1, 1, 1])
title('sxx')

