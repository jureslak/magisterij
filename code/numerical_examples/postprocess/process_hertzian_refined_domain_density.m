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

% testidx = find(x.^2 + y.^2 < (1000*b)^2);% & y.^2 > (b/2)^2);
% 
% errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
% erryy = max(max(abs(syy(testidx) - asyy(testidx))));
% errxy = max(max(abs(sxy(testidx) - asxy(testidx))));
% 
% err = max([errxx, erryy, errxy])/p0

lim = 5*b;
interesting = -lim < x & x < lim & -lim < y;

xi = double(x(interesting));
yi = double(y(interesting));
sxxi = sxx(interesting);
syyi = syy(interesting);
sxyi = sxy(interesting);
vmi = von_mises(sxxi, syyi, sxyi);
[asxxi, asyyi, asxyi] = hertzian_analytical(xi, yi, b, p0);
avmi = von_mises(asxxi, asyyi, asxyi);

NS = KDTreeSearcher([x', y']);
[idx, d] = knnsearch(NS, [x', y'], 'Distance', 'euclidean', 'K', 2);
avgd = double(d(:, 2:end));
avgd = avgd(interesting) / max(avgd);

save([plotdatapath 'hertzian_refined_domain_density.mat'], 'b', 'p0',...
     'xi', 'yi', 'avgd');