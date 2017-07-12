prepare
datafile = [datapath 'hertzian_matrix_example.h5'];
M = spconvert(h5read(datafile, '/M')');
save([plotdatapath 'hertzian_matrix_example.mat'], 'M')
cspy(M)