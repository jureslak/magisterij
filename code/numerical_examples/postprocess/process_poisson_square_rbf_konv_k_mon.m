prepare
datafile = [datapath 'poisson_square_implicit_rbf_konv_k_mon.h5'];
info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
sigmas = zeros(typenum, 1);

for i = 1:simnum
    pos = h5read(datafile, [info.Groups(1).Groups(i).Name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    anal = poisson_square_analytical(x, y)';
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        sol = double(h5read(datafile, [name '/sol']));
        N = h5readatt(datafile, name, 'N');
        time = h5readatt(datafile, name, 'timetotal');

        err = max(max(abs(sol - anal)));
        cutoff = h5read(datafile, [name, '/cutoff']);
        cutoff = mean(reshape(cutoff, [2, length(cutoff)/2]));
        
        data(j, i, 1) = N;
        data(j, i, 2) = err;
        data(j, i, 3) = mean(cutoff);
        
        sigmas(j) = str2double(name(9:11));
    end
    fprintf('point %d/%d \r', i, simnum);
end

names = cell(typenum, 1);
for i = 1:typenum
    names{i} = info.Groups(i).Name;
end

save([plotdatapath 'poisson_square_rbf_konv_k_mon.mat'], 'data', 'names', 'typenum', 'simnum', 'sigmas')