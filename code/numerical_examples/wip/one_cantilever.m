prepare
datafile = [datapath 'cantilever_convergence2.h5'];

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
   
[asxx, asyy, asxy, au, av] = cantilever_beam_analytical(x, y, P, L, D, E, nu);
asv = von_mises(asxx, asyy, asxy);

% M = spconvert(h5read(datafile, [name '/matrix'])');
% rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% max(max(displ2))

%%
close all

sf1 = setfig('b1');
% scatter(x, y, 5, sxx, 'filled')
f = 1e5;
scontour(x+f*u, y+f*v, sxx, 100, 100, 20);
colorbar
daspect([1, 1, 1])
title('sxx')

sf2 = setfig('b2');
% scatter(x, y, 5, syy, 'filled')
scontour(x, y, syy, 100, 100, 20);
colorbar
daspect([1, 1, 1])
title('syy')

sf3 = setfig('b3');
% scatter(x, y, 5, sxy, 'filled')
scontour(x, y, sxy, 100, 100, 20);
colorbar
title('sxy')
daspect([1, 1, 1])

errxx = max(max(abs(sxx - asxx)));
erryy = max(max(abs(syy - asyy)));
errxy = max(max(abs(sxy - asxy)));
errstress =  max([errxx, erryy, errxy])

Mu = max(sqrt(au.^2 + av.^2));
erru = max(max(abs(u - au)))/Mu;
errv = max(max(abs(v - av)))/Mu;
erruv = max([erru, errv])

f2 = setfig('bo11');
top = find(y == max(y));
[~, I] = sort(x(top));
top = top(I);

plot(x(top), u(top), 'o');
plot(x(top), au(top), '-');

plot(x(top), v(top), 'o');
plot(x(top), av(top), '-');

% plot(x(top), syy(top), 'o');
% plot(x(top), asyy(top), '-');

% plot(x(top), sxy(top), 'o');
% plot(x(top), asxy(top), '-');


% plot(x(top), sv(top), 'o');
% plot(x(top), asv(top), '-');

% legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
%        'von mises analytical', 'Location', 'SE')
legend('u', 'au', 'v', 'av', 'syy', 'asyy', 'sxy', 'asxy', 'Location', 'SE')
% xlim([-3, 3])


f3 = setfig('bo15');
[mid, xmid] = closest(pos, 1, L);
plot(u(mid), y(mid), 'o');
plot(au(mid), y(mid), '-');

plot(v(mid), y(mid), 'o');
plot(av(mid), y(mid), '-');

% plot(sxy(mid), y(mid), 'o');
% plot(asxy(mid), y(mid), '-');

% plot(syy(mid), y(mid), 'o');
% plot(syy(mid), y(mid), '-');

legend('u', 'au', 'v', 'av', 'sxy', 'asxy', 'syy', 'asyy', 'Location', 'SE')

f4 = setfig('bo16');
top = find(y == min(y));
[~, I] = sort(x(top));
top = top(I);

plot(x(top), u(top), 'o');
plot(x(top), au(top), '-');

plot(x(top), v(top), 'o');
plot(x(top), av(top), '-');

% plot(x(top), syy(top), 'o');
% plot(x(top), asyy(top), '-');

% plot(x(top), sxy(top), 'o');
% plot(x(top), asxy(top), '-');

% legend('u', 'au', 'v', 'av', 'Location', 'SE')
% 
f5 = setfig('bo12');
[mid, xmid] = closest(pos, 1, 0);
plot(u(mid), y(mid), 'o');
plot(au(mid), y(mid), '-');

plot(v(mid), y(mid), 'o');
plot(av(mid), y(mid), '-');

legend('u', 'au', 'v', 'av', 'Location', 'SE')

% plot(x(top), sv(top), 'o');
% plot(x(top), asv(top), '-');

% legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
%        'von mises analytical', 'Location', 'SE')
% legend('u', 'au', 'v', 'av', 'Location', 'SE')

% ylim([-3, 3])