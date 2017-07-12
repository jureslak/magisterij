prepare
load([plotdatapath 'poisson_square_implicit.mat'])

keys = {'/gau13','/gau5','/gau9','/imq13','/imq5','/imq9','/mon5','/mon9','/mq13','/mq5','/mq9', '/mon6'};
vals = {'G13','G5','G9','IMQ13','IMQ5','IMQ9', 'MON5', 'MON9', 'MQ13','MQ5','MQ9', 'MON6'};
namemap = containers.Map(keys, vals);

legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = namemap(names{i});
end

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
colors = {
    [1,123,206]/256,
    [1,56,147]/256,
%     [0,98,199]/256,
    [253,94,91]/256,
    [137,1,1]/256,
%     [238,16,31]/256,
    [255,207,0]/256,
    [255,169,0]/256,
    [223,129,9]/256,
    [57,172,55]/256,
    [15,85,48]/256,
%     [19,131,49]/256
};

Ns = data(1, :, 1);

best = polyfit(log(Ns), log(data(5, :, 2)), 1);
f1 = setfig('b1');
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
text(0.4, 0.43, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([40, 10^5])
legend(legendvals)

f2 = setfig('b2', 'Visible', 'on');
hold off
hold on
for i = 1:typenum
    plot(Ns, data(i, :, 3), [markers{i}, '-'], 'Color', colors{i})
end
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([40, 10^5])
legend(legendvals, 'Location', 'NW')
xlabel('$N$')
ylabel('\v{c}as [$s$]')

exportfig(f1, [imagepath 'poisson_square_convergence'], '-pdf')
exportfig(f2, [imagepath 'poisson_square_time'], '-pdf')