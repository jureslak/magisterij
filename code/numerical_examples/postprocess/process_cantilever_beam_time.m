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
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        data(j, i, 1) = N;
        
        time(i, :) = [h5readatt(datafile, name, 'time_domain');
                      h5readatt(datafile, name, 'time_shapes');
                      h5readatt(datafile, name, 'time_construct');
                      h5readatt(datafile, name, 'time_lut');
                      h5readatt(datafile, name, 'time_solve');
                      h5readatt(datafile, name, 'time_post')];
    end
    fprintf('point %d/%d %s\r', i, simnum, name);
end

save([plotdatapath 'cantilever_beam_time.mat'], 'data', 'time');