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
times = errs;
timespart = errs;
timem = errs;
Ns = zeros(simnum, 1);
for i = 1:simnum
    grp = info.Groups(i);
    name = grp.Name;
    
    xx = h5read(datafile, [name '/pos']);
    yy = h5read(datafile, [name '/sol']);
    yyfdm = h5read(datafile, [name '/solfdm']);
    time = h5readatt(datafile, name, 'time');
    timefdm = h5readatt(datafile, name, 'timefdm');
    timepart = h5readatt(datafile, name, 'timepart');

    timem(i) = timefdm;
    times(i) = time;
    timespart(i) = timepart;
    
    errs(i) = max(abs(correct(xx)' - yy));
    errsfdm(i) = max(abs(correct(linspace(0, 1, length(xx)))' - yyfdm));
    Ns(i) = length(xx);

    %if i > 50, break, end
end

ok = find(Ns < 1e3);
best = polyfit(log(Ns(ok)), log(errs(ok)), 1);

f1 = setfig('b1');
plot(Ns, errs, 'o-')
plot(Ns, errsfdm, '^-')
plot(Ns, exp(best(2))*Ns.^best(1), '--b');
text(0.4, 0.52, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([0, 2*10^5])
ylim([1e-11 1e0])
set(gca, 'XScale', 'log', 'YScale', 'log')
legend('MLSM', 'FDM', 'trend')

f2 = setfig('b2');
plot(Ns, times, 'o-')
plot(Ns, times - timespart, 'x-')
plot(Ns, timem, '^-');
ylim([0 0.5])
xlabel('$N$', 'FontSize', 14)
ylabel('\v{c}as [$s$]', 'FontSize', 14)
%set(gca, 'XScale', 'log', 'YScale', 'log')
yyaxis right
ylabel('razmerje')
set(gca,'Position', [.1 .1 .83 .85])
plot(Ns, times ./ timem, 's-');
legend('MLSM-skupaj', 'MLSM-re\v{s}evanje', 'FDM', 'razmerje', 'Location', 'E')

figure(f1);
figure(f2);

% save

exportfig(f1, '../../../images/lap1d_convergence', '-pdf')
exportfig(f2, '../../../images/lap1d_times', '-pdf')