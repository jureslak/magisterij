prepare 

datafile = [datapath 'hertzian_domain_too_small.h5'];
info = h5info(datafile);

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');
max_height = h5readatt(datafile, '/', 'max_height');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);
dxs = zeros(simnum, 1);

tops = cell(simnum, 2);

% f2 = setfig('b1');
% f3 = setfig('b3');
% colors = colormap;
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;
    
    for j = 1:typenum
%         if j > 7, break, end
        grp = info.Groups(j).Groups(i);
        name = grp.Name
        
        N = h5readatt(datafile, name, 'N');
        nglob = h5readatt(datafile, name, 'nglob');

        height = h5readatt(datafile, name, 'height');
        pos = h5read(datafile, [name '/pos']);
        x = pos(1, :);
        y = pos(2, :);

        testidx = find(x.^2 + y.^2 < (10*b)^2);
        [asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);
        asv = von_mises(asxx, asyy, asxy);
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);
        sv = von_mises(sxx, syy, sxy);

        errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
        erryy = max(max(abs(syy(testidx) - asyy(testidx))));
        errxy = max(max(abs(sxy(testidx) - asxy(testidx))));

        data(j, i, 1) = N * (max_height*b / height)^2;
        data(j, i, 2) = max([errxx, erryy, errxy])/p0;
        data(j, i, 3) = height / b;        
    end

    dxs(i) = nglob/max_height;
    
    fprintf('point %d/%d %s\r', i, simnum, name);
end

%%
close all
markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
% legendvals = {'MON9'};
colors = {
    [1,123,206]/256,
    [255,169,0]/256,
    [253,94,91]/256,
    [19,131,49]/256,
    [255,169,0]/256,
    [15,85,48]/256+0.1,
    [57,172,55]/256,
    [1,123,206]/256,
    [255,207,0]/256,
    [223,129,9]/256,
    [19,131,49]/256
};
names = cell(simnum, 1);
for i = 1:simnum
    names{i} = sprintf('%.2f nodes / $b$', dxs(i));
end

f1 = setfig('b1');
for i = 1:simnum
    plot(data(:, i, 3), data(:, i, 2), [markers{mod(i, 7)+1}, '-'], 'Color', colors{i})
end
% set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$H/b$')
ylabel('$L_\infty$ napaka')
% xlim([1e2, 2e6])
ylim([0.05, 0.18])
% set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(names)