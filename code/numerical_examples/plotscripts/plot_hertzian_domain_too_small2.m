prepare 
load([plotdatapath 'hertzian_domain_too_small2.mat']);

[typenum, simnum, ~] = size(data);

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
% legendvals = {'MON9'};
colors = {
    [1,123,206]/256,
    [255,169,0]/256,
    [253,94,91]/256,
    [19,131,49]/256,
    [255,169,0]/256,
    [15,85,48]/256+0.1,
    [57,172,55]/256,
    [1,123,206]/256,
    [255,207,0]/256,
    [223,129,9]/256,
    [19,131,49]/256
};
names = cell(simnum, 1);
for i = 1:simnum
    names{i} = sprintf('%.2f to\\v{c}k / $b$', dxs(i));
end

f1 = setfig('b1');
for i = 1:simnum
    plot(data(:, i, 3), data(:, i, 2), [markers{mod(i, 7)+1}, '-'], 'Color', colors{i})
end
% set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$H/b$')
ylabel('$L_\infty$ napaka')
% xlim([1e2, 2e6])
ylim([0.05, 0.18])
% set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(names)

exportfig(f1, [imagepath 'hertzian_domain_too_small2'], '-pdf')