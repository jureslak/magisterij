close all
clear
datafile = '../data/poisson_neumann.h5';
info = h5info(datafile);
simnum = length(info.Groups);

evalpt = 500;
[X,Y] = meshgrid(linspace(0, 1, evalpt), linspace(0, 1, evalpt));

name = info.Groups(1).Name;
pos = double(h5read(datafile, [name '/pos']));
x = pos(1, :);
y = pos(2, :);
sol = double(h5read(datafile, [name '/sol']));
Zp = griddata(x, y, sol, X, Y);

errs = zeros(simnum, 1);
Ns = zeros(simnum, 1);
for i = 2:simnum
    grp = info.Groups(i);
    name = grp.Name;
    
    pos = double(h5read(datafile, [name '/pos']));
    sol = double(h5read(datafile, [name '/sol']));
    
    x = pos(1, :);
    y = pos(2, :);
    Z = griddata(x, y, sol, X, Y);
    
    
    Ns(i) = length(pos);
    errs(i) = max(max(abs(Z-Zp)));
    
    Zp = Z;
    i=i
end

%%

close all
f1 = setfig('b1');
plot(Ns, errs, 'o-')
xlabel('$N$')
ylabel('$L_\infty$ razlika med zaporednima re\v{s}itvama')
%xlim([0, 2*10^5])
%ylim([1e-11 1e0])
set(gca, 'XScale', 'log', 'YScale', 'log')
figure(f1);

