close all
clear
datafile = '../data/convergence_test_neumann.h5';
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

ok = find(Ns < 1e3);
best = polyfit(log(Ns(ok)), log(errs(ok)), 1);
ok = find(Ns < 30);
best2 = polyfit(log(Ns(ok)), log(errs5(ok)), 1);


f1 = setfig('b1');
plot(Ns, errs, 'o-')
plot(Ns, errsfdm, '^-')
plot(Ns, errs_float, 's-')
plot(Ns, errsfdm_float, 'p-')
plot(Ns, errs5, '*-')
plot(Ns, exp(best(2))*Ns.^best(1), '--k');
plot(Ns, exp(best2(2))*Ns.^best2(1), '--k');

text(0.136, 0.45, sprintf('$k = %.2f$', best2(1)), 'Units', 'normalized')
text(0.4, 0.45, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([0, 10^5])
ylim([1e-9 1e-1])
set(gca, 'XScale', 'log', 'YScale', 'log')
legend('MLSM double', 'FDM double', 'MLSM float', 'FDM float', 'MLSM double vi\v{s}ji red')


f2 = setfig('b2');
plot(Ns, 1000*times, 'o-')
plot(Ns, 1000*timesfdm, '^-');
plot(Ns, time_float, '^-');
plot(Ns, timefdm_float, '^-');
xlim([0, 1e5])
plot(Ns, 1000*time5, 'p-');
plot(Ns, 1000*times - 1000*timespart, 'x-')
ylim([1e-5 1e4])
xlabel('$N$', 'FontSize', 14)
ylabel('\v{c}as [ms]', 'FontSize', 14)
set(gca, 'XScale', 'log', 'YScale', 'log')
yyaxis right
ylabel('razmerje')
set(gca,'Position', [.1 .1 .83 .85])
grey = 0.5*ones(3, 1);
set(gca, 'YColor', grey)
plot(Ns, times ./ timesfdm, 's-.', 'Color', grey);
plot(Ns, time5 ./ times, '*-', 'Color', grey);
legend('MLSM double -- skupaj', 'FDM double', 'MLSM float', 'FDM float',...
       'MLSM double vi\v{s}ji red -- skupaj',...
       'MLSM double -- re\v{s}evanje', 'razmerje MLSM / FDM',...
       'razmerje MLSM vi\v{s}ji red / MLSM', 'Location', 'NW')
figure(f1);
figure(f2);

% save

exportfig(f1, '../../../images/lap1d_convergence', '-pdf')
exportfig(f2, '../../../images/lap1d_times', '-pdf')
