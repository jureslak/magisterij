clear
datafile = '../data/poisson_square_implicit_rbf_konv_k_mon.h5';
info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
sigmas = zeros(typenum, 1);

for i = 1:simnum
    pos = h5read(datafile, [info.Groups(1).Groups(i).Name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    anal = poisson_square_analytical(x, y)';
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        sol = double(h5read(datafile, [name '/sol']));
        N = h5readatt(datafile, name, 'N');
        time = h5readatt(datafile, name, 'timetotal');

        err = max(max(abs(sol - anal)));
        cutoff = h5read(datafile, [name, '/cutoff']);
        cutoff = mean(reshape(cutoff, [2, length(cutoff)/2]));
        
        data(j, i, 1) = N;
        data(j, i, 2) = err;
        data(j, i, 3) = mean(cutoff);
        
        sigmas(j) = str2double(name(9:11));
    end
    fprintf('point %d/%d \r', i, simnum);
end

%%

close all
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

exportfig(f1, '../../../images/poisson_square_rbf_konv_k_mon', '-pdf')
% exportfig(f2, '../../../images/poisson_square_time', '-pdf')