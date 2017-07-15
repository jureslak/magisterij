prepare
load([plotdatapath, 'hertzian_refine_levels_convergence.mat'])

[typenum, simnum, ~] = size(data);

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
lines = {'-', '--', '-', '--', '-', '--', '-'};
legendvals = {'nivo 0', 'nivo 1', 'nivo 2', 'nivo 3', 'nivo 4', 'nivo 5', 'nivo 6'};
colors = {
    [1,123,206]/256, % blue
    [255,169,0]/256, % yellow
    [253,94,91]/256, % red
    [57,172,55]/256, % green
    [115,0,171]/256, % purple
    [100,100,100]/256, % gray
    [115,64,33]/256, % brown
    [223,129,9]/256,
};


f1 = setfig('b1');
for i = 1:typenum
    plot(data(i, :, 1), data(i, :, 2), [markers{i}, lines{i}], 'Color', colors{i})
end

% best fit
fitidx = 3;
fitrange = 9:22;
best = polyfit(log(data(fitidx, fitrange, 1)), log(data(fitidx, fitrange, 2)), 1);
plot(data(fitidx, fitrange, 1), exp(best(2))*data(fitidx, fitrange, 1).^best(1), ':k', 'LineWidth', 1.5);
text(0.853, 0.414, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')

set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([8e3, 2e6])
ylim([-inf, 5e-2])
set(gca, 'YTick', [5e-3, 1e-2, 5e-2])
legendvals{end+1} = 'trend';
legend(legendvals)

exportfig(f1, [imagepath 'hertzian_refine_levels_convergence'], '-pdf')