prepare
datafile = [datapath 'poisson_weird3d.h5'];
pos = h5read(datafile, '/calc/pos');
sol = h5read(datafile, '/calc/sol');
save([plotdatapath 'poisson_weird.mat'], 'pos', 'sol')