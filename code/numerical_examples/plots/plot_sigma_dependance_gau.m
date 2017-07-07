filename = '../data/poisson_square_implicit_sigma_scan.h5';
info = h5info(filename);
groups = info.Groups;
sbn = length(groups);
swn = length(groups(1).Groups);

sbs = zeros(sbn, 1);
sws = zeros(swn, 1);

pos = h5read(filename, sprintf('%s/pos', groups(1).Groups(1).Name));
N = length(pos);
x = pos(1, :);
y = pos(2, :);
anal = poisson_square_analytical(x, y);

cutoff = zeros(sbn, swn);
err = zeros(sbn, swn);
for i = 1:sbn
    for j = 1:swn
        name = groups(i).Groups(j).Name;
        fprintf('%s\n', name)
        sbs(i) = h5readatt(filename, name, 'sigmaB');
        sws(j) = h5readatt(filename, name, 'sigmaW');
        cutoff(i, j) = h5readatt(filename, name, 'cutoff');

        sol = h5read(filename, sprintf('%s/sol', name))';
        err(i, j) = max(abs(sol - anal));
    end
end

%%
close all
xlables = cell(length(sbs), 1);
for i = 1:length(sbs)
    if mod(i, 2) == 0
        xlabels{i} = sprintf('%.1f', sbs(i));
    else
        xlabels{i} = '';
    end
end   
ylables = cell(length(sws), 1);
for i = 1:length(sws)
    if mod(i, 2) == 1 || sws(i) > 1
        ylabels{i} = sprintf('%.2f', sws(i));
    else
        ylabels{i} = '';
    end
end

f1 = setfig('b1', 'Visible', 'on');
hold off
% set(gca,'Position', [.1 .15 .85 .80])
[W,B] = meshgrid(sws, sbs);
err(err > 1e-1) = 1e-1;
contour(B, W, log10(err), 'Fill', 'on')
title('$\log_{10}(L_\infty$ napaka)')
c = colorbar;
xlabel(c, sprintf('$n = m = 9, N = %d$, Gaussove bazne funkcije', N), 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\sigma_b / r_\chi$')
set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
ylabel('$\sigma_w / r_\chi$')
set(gca,'YTick', sws,'YTickLabel', ylabels)
fmid = gcf;
% set(gca, 'XScale', 'log', 'YScale', 'log')


f2 = setfig('b2', 'Visible', 'on');
hold off
contour(B, W, cutoff, 'Fill', 'on')
title('\v{S}tevilo odrezanih singularnih vrednosti')
c = colorbar;
xlabel(c, sprintf('$n = m = 9, N = %d$, Gaussove bazne funkcije', N), 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\sigma_b / r_\chi$')
set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
ylabel('$\sigma_w / r_\chi$')
set(gca,'YTick', sws,'YTickLabel', ylabels)
fcut = gcf;

figure(f1)

exportfig(f1, '../../../images/poisson_square_sigma_depedence_error', '-png');
exportfig(f2, '../../../images/poisson_square_sigma_depedence_cutoff', '-png');