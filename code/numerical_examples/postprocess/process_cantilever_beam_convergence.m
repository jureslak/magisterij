prepare
datafile = [datapath 'cantilever_convergence_jarjar.h5'];
info = h5info(datafile);

P = h5readatt(datafile, '/', 'P');
D = h5readatt(datafile, '/', 'D');
L = h5readatt(datafile, '/', 'L');
I = h5readatt(datafile, '/', 'I');
E = h5readatt(datafile, '/', 'E');
nu = h5readatt(datafile, '/', 'v');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');
    pos = h5read(datafile, [name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    
    [asxx, asyy, asxy, au, av] = cantilever_beam_analytical(x, y, P, L, D, E, nu);
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);
        
        displ  = h5read(datafile, [name '/disp']);
        u = displ(1, :);
        v = displ(2, :);
        
        erru = max([max(abs(u - au)), max(abs(v - av))]);
        Mu = max([max(abs(au)), max(abs(av))]);

        errxx = max(max(abs(sxx - asxx)));
        erryy = max(max(abs(syy - asyy)));
        errxy = max(max(abs(sxy - asxy)));
        Ms = max([max(abs(sxx)) max(abs(syy)) max(abs(sxy))]);

        data(j, i, 1) = N;
        data(j, i, 2) = erru / Mu;
        data(j, i, 4) = max([errxx, erryy, errxy]) / Ms;
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');
        
        time(i, :) = [h5readatt(datafile, name, 'time_domain');
                      h5readatt(datafile, name, 'time_shapes');
                      h5readatt(datafile, name, 'time_construct');
                      h5readatt(datafile, name, 'time_lut');
                      h5readatt(datafile, name, 'time_solve');
                      h5readatt(datafile, name, 'time_post')];
    end
    fprintf('point %d/%d %s\r', i, simnum, name);
end

save([plotdatapath 'cantilever_beam_convergence.mat'], 'data');