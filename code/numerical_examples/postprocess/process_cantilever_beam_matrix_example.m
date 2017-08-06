prepare
datafile = [datapath 'cantilever_beam_matrix_example.h5'];
M = spconvert(h5read(datafile, '/M')');
save([plotdatapath 'cantilever_beam_matrix_example.mat'], 'M')