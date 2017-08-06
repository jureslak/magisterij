prepare
load([plotdatapath 'cantilever_beam_convergence.mat'])
[typenum, simnum, ~] = size(data);
markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
legendvals = {'G9 -- pomik ($e_\infty(\vec{u})$)', 'G9 -- ($e_\infty(\sigma)$)', 'MON9 -- pomik ($e_\infty(\vec{u})$)',...
    'MON9 -- ($e_\infty(\sigma)$)'};
colors = {
    [1,123,206]/256,
    [255,169,0]/256,
%     [1,56,147]/256,
%     [0,98,199]/256,
    [253,94,91]/256,
    [255,169,0]/256,
%     [137,1,1]/256,
%     [238,16,31]/256,
    [15,85,48]/256+0.1,
    [57,172,55]/256,
    [1,123,206]/256,
    [255,207,0]/256,
    [223,129,9]/256,
    [19,131,49]/256
};

f1 = setfig('b1');
Ns = data(1, :, 1);
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{2*i}, '-'], 'Color', colors{i})
    plot(Ns, data(i, :, 4), [markers{2*i+1}, '--'], 'Color', colors{i})
end

fitrange = 1:40;
fitx = Ns(fitrange);
best = polyfit(log(fitx), log(data(2, fitrange, 2)), 1);
% plot(fitx, exp(best(2))*fitx.^best(1), ':k');
text(0.775, 0.131, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')

best2 = polyfit(log(fitx), log(data(2, fitrange, 4)), 1);
text(0.873, 0.29, sprintf('$k = %.2f$', best2(1)), 'Units', 'normalized')

set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([1e2, 1.1e6])
ylim([5e-6, 1e0])
% set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(legendvals)

exportfig(f1, [imagepath 'cantilever_beam_convergence'], '-pdf')