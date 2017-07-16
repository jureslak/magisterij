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

save([plotdatapath 'hertzian_domain_too_small2.mat'], 'data', 'dxs');
