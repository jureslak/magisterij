prepare
datafile = [datapath 'hertzian_domain_too_small.h5'];

name = '/mon_h00005.000000/calc0014_0113';

pos = h5read(datafile, [name '/pos']);
stress = h5read(datafile, [name '/stress']);
N = h5readatt(datafile, name, 'N');

x = pos(1, :);
y = pos(2, :);

sxx = stress(1, :);
syy = stress(2, :);
sxy = stress(3, :);
sv = von_mises(sxx, syy, sxy);

displ  = h5read(datafile, [name '/disp']);
u = displ(1, :);
v = displ(2, :);
    
max(max(displ))

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
[asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);
asv = von_mises(asxx, asyy, asxy);

% M = spconvert(h5read(datafile, [name '/matrix'])');
% rhs = h5read(datafile, [name '/rhs']);

% displ2 = M \ rhs;
% max(max(displ2))

%%
close all

testidx = find(x.^2 + y.^2 < (10*b)^2); % & y.^2 > b^2);
% testidx = find(x.^2 + y.^2 < (10*b)^2 & y.^2 > (0.1*b)^2);
%testidx = find(x.^2 + y.^2 < (1000*b)^2);

f1 = setfig('b1');
scatter(x(testidx)/b, y(testidx)/b, 5, log10(abs(sxx(testidx) - asxx(testidx))/p0), 'filled')
colorbar
daspect([1, 1, 1])


errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
erryy = max(max(abs(syy(testidx) - asyy(testidx))));
errxy = max(max(abs(sxy(testidx) - asxy(testidx))));
err =  max([errxx, erryy, errxy])/p0
errvm = max(abs(sv(testidx) - asv(testidx)))/p0

f2 = setfig('b2');
top = find(y == max(y));
[~, I] = sort(x(top));
top = top(I);

plot(x(top)/b, sxx(top)/p0, 'o');
plot(x(top)/b, asxx(top)/p0, '-');

plot(x(top)/b, syy(top)/p0, 'o');
plot(x(top)/b, asyy(top)/p0, '-');

plot(x(top)/b, sxy(top)/p0, 'o');
plot(x(top)/b, asxy(top)/p0, '-');

plot(x(top)/b, sv(top)/p0, 'o');
plot(x(top)/b, asv(top)/p0, '-');

legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
       'von mises analytical', 'Location', 'SE')
xlim([-3, 3])


f4 = setfig('b4');
[mid, xmid] = closest(pos, 1, 0);
plot(sxx(mid)/p0, y(mid)/b, 'o');
plot(asxx(mid)/p0, y(mid)/b, '-');

plot(syy(mid)/p0, y(mid)/b, 'o');
plot(asyy(mid)/p0, y(mid)/b, '-');

plot(sxy(mid)/p0, y(mid)/b, 'o');
plot(asxy(mid)/p0, y(mid)/b, '-');

plot(sv(mid)/p0, y(mid)/b, 'o');
plot(asv(mid)/p0, y(mid)/b, '-');

legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
       'von mises analytical', 'Location', 'SE')

% ylim([-3, 3])