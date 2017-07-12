load([plotdatapath 'poisson_sigma_dependence_gau.mat'])

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

exportfig(f1, [imagepath 'poisson_square_sigma_depedence_error'], '-r300', '-png', '-opengl');
exportfig(f2, [imagepath 'poisson_square_sigma_depedence_cutoff'], '-r300', '-png', '-opengl');