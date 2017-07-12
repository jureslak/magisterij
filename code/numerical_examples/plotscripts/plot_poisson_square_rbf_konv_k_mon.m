prepare
load([plotdatapath 'poisson_square_rbf_konv_k_mon.mat'])

legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = sprintf('$%.0f\\,r_\\chi$', sigmas(i));
end
legendvals{end} = '$\infty$';

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
styles = {'-', '-.'};
colors = {
    [1,56,147]/256,
    [1,123,206]/256,
    [0,98,199]/256,
    [137,1,1]/256,
    [253,94,91]/256,
    [238,16,31]/256,
    [255,207,0]/256,
    [255,169,0]/256,
    [223,129,9]/256,
    [15,85,48]/256,
    [57,172,55]/256,
    [19,131,49]/256
};
setfig('a1', 'Visible', 'off');
cmap = colormap('parula');
cmap = cmap(1:5:length(cmap), :);
colors = mat2cell(cmap, ones(length(cmap), 1), 3);

Ns = data(1, :, 1);

f1 = setfig('b1');
for i = 1:typenum
    last = plot(Ns, data(i, :, 2), ...
        [markers{mod(i, length(markers))+1},...
         styles{mod(i, length(styles))+1}], 'Color', colors{i});
end
set(last, 'Color', [0.3 0.3 0.3], 'LineStyle', '-.', 'Marker', 'x');
% uistack(last, 'bottom')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([40, 10^5])
ylim([10e-8, 1e-2])
hleg = legend(legendvals, 'Location', 'SW');
title(hleg, '$\sigma_B$');

exportfig(f1, [imagepath 'poisson_square_rbf_konv_k_mon'], '-pdf')
% exportfig(f2, '../../../images/poisson_square_time', '-pdf')