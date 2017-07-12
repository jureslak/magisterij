prepare
load([plotdatapath 'lap1d.mat'])

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
plot(Ns, 1000*time_float, '^-');
plot(Ns, 1000*timefdm_float, '^-');
xlim([0, 1e5])
plot(Ns, 1000*time5, 'p-');
plot(Ns, 1000*times - 1000*timespart, 'x-')
ylim([1e-5 1e5])
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

exportfig(f1, [imagepath 'lap1d_convergence'], '-pdf')
exportfig(f2, [imagepath 'lap1d_times'], '-pdf')