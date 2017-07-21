prepare
datafile = [datapath 'poisson_square_parallel.h5'];

info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 2);
time = zeros(typenum, simnum, 5);
threads = zeros(typenum, 1);
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');

    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        data(j, i, 1) = N;
        data(j, i, 2) = h5readatt(datafile, name, 'time_total');
        threads(j) = h5readatt(datafile, name, 'thread_num');
        
        time(j, i, :) = [h5readatt(datafile, name, 'time_domain');
                         h5readatt(datafile, name, 'time_shapes');
                         h5readatt(datafile, name, 'time_construct')
                         h5readatt(datafile, name, 'time_compute')
                         h5readatt(datafile, name, 'time_solve');];
    end

    fprintf('point %d/%d \r', i, simnum);
end


save([plotdatapath 'poisson_square_parallel.mat'], 'data', 'time', 'threads');