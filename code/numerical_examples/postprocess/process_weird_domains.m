prepare

datafile = [datapath 'unusual_domains.h5'];

lpos = hdf5read(datafile, '/luna/pos');
ltype = hdf5read(datafile, '/luna/types');
dpos = hdf5read(datafile, '/drevo/pos');
dtype = hdf5read(datafile, '/drevo/types');
cpos = hdf5read(datafile, '/cev/pos');
ctype = hdf5read(datafile, '/cev/types');

save([plotdatapath 'weird_domains.mat'], 'lpos', 'ltype',...
      'cpos', 'ctype', 'dpos', 'dtype');
