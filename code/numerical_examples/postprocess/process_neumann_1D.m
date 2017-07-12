prepare
datafile = [datapath 'convergence_test_neumann.h5'];
info = h5info(datafile);
simnum = length(info.Groups);

correct = @(x) cos(1)*x - sin(x);
xanal = linspace(0, 1, 1000);
yanal = correct(xanal);

errs = zeros(simnum, 1);
errsfdm = errs;
errs_float = errs;
errsfdm_float = errs;
errs5 = errs;
times = errs;
timespart = errs;
timesfdm = errs;
timefdm_float = errs;
time_float = errs;
time5 = errs;
Ns = zeros(simnum, 1);
for i = 1:simnum
    grp = info.Groups(i);
    name = grp.Name;

    xx = h5read(datafile, [name '/pos']);
    if (length(xx) < 500)
        yy_float = h5read(datafile, [name '/sol_float']);
        yyfdm_float = h5read(datafile, [name '/solfdm_float']);
        timefdm_float(i) = h5readatt(datafile, name, 'timefdm_float');
        time_float(i) = h5readatt(datafile, name, 'time_float');
        errs_float(i) = max(abs(correct(xx)' - yy_float));
        errsfdm_float(i) = max(abs(correct(linspace(0, 1, length(xx)))' - yyfdm_float));
    else
        time_float(i) = nan;
        timefdm_float(i) = nan;
    end

    yy = h5read(datafile, [name '/sol']);
    yyfdm = h5read(datafile, [name '/solfdm']);
    yy5 = h5read(datafile, [name '/sol_s5']);

    timesfdm(i) = h5readatt(datafile, name, 'timefdm');
    times(i) = h5readatt(datafile, name, 'time');
    time5(i) = h5readatt(datafile, name, 'time_s5');
    timespart(i) =  h5readatt(datafile, name, 'timepart');

    errs(i) = max(abs(correct(xx)' - yy));
    errsfdm(i) = max(abs(correct(linspace(0, 1, length(xx)))' - yyfdm));
    errs5(i) = max(abs(correct(xx)' - yy5));

    Ns(i) = length(xx);

    %if i > 50, break, end
end

save([plotdatapath 'lap1d.mat'], 'errs', 'errsfdm', 'errs_float', 'errs5',...
     'errsfdm_float', 'Ns', 'times', 'timesfdm', 'timefdm_float', 'time_float',...
     'timespart', 'time5');