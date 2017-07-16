prepare
datafile = [datapath 'hertzian_refined_convergence_jarjar.h5'];
info = h5info(datafile);

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

% f2 = setfig('b1');
% f3 = setfig('b3');
% colors = colormap;
for i = 1:simnum
    if i > 22, break, end
    if i < 9, continue, end
   
    for j = 1:typenum
        if i == 22 && j > 5, continue, end
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        N = h5readatt(datafile, name, 'N');
        pos = h5read(datafile, [name '/pos']);
        x = pos(1, :);
        y = pos(2, :);

        testidx = find(x.^2 + y.^2 < (1000*b)^2);% & y.^2 > (b/2)^2);
        [asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);

        errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
        erryy = max(max(abs(syy(testidx) - asyy(testidx))));
        errxy = max(max(abs(sxy(testidx) - asxy(testidx))));

        data(j, i, 1) = N;
        data(j, i, 2) = max([errxx, erryy, errxy])/p0;
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');
        
        time(i, :) = [h5readatt(datafile, name, 'time_domain');
                      h5readatt(datafile, name, 'time_shapes');
                      h5readatt(datafile, name, 'time_construct');
                      h5readatt(datafile, name, 'time_lut');
                      h5readatt(datafile, name, 'time_solve');
                      h5readatt(datafile, name, 'time_post')];
    end
    fprintf('point %d/%d \r', i, simnum);
end

save([plotdatapath, 'hertzian_refine_levels_convergence.mat'],...
     'data');
