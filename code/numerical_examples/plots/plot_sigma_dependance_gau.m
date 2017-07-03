close all
sbs = dlmread('../plotdata/sigmaB_gau.txt');
sws = dlmread('../plotdata/sigmaW_gau.txt');
cutoff = dlmread('../plotdata/number_cuttoff_svd_gau.txt');
% convergence = dlmread('plotdata/bicgstab_error_estimate_gau.txt');
% iters = dlmread('plotdata/bicgstab_iterations_gau.txt');
err = dlmread('../plotdata/error_gau.txt');

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
    if mod(i, 1) == 0
        ylabels{i} = sprintf('%.2f', sws(i));
    else
        ylabels{i} = '';
    end
end

setfig('b1', 'Visible', 'on')
hold off
% set(gca,'Position', [.1 .15 .85 .80])
[W,B] = meshgrid(sws, sbs);
err(err > 10) = 10;
contour(B, W, log10(err), 'Fill', 'on')
title('$\log_{10}(L_\infty$ napaka), Gaussove bazne funkcije.')
c = colorbar;
xlabel(c, '$n = m = 9, N = 3751, \Delta x = \Delta y = 3.33\cdot 10^{-4}$', 'interpreter', 'latex', 'fontsize', 12)
xlabel('$\sigma_b / r_\chi$')
set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
ylabel('$\sigma_w / r_\chi$')
set(gca,'YTick', sws,'YTickLabel', ylabels)
fmid = gcf;

setfig('b2', 'Visible', 'on')
hold off
contour(B, W, cutoff, 'Fill', 'on')
title('\v{S}tevilo odrezanih singularnih vrednosti, Gaussove bazne funkcije.')
c = colorbar;
xlabel(c, '$n = m = 9, N = 3751, \Delta x = \Delta y = 3.33\cdot 10^{-4}$')
xlabel('$\sigma_b / r_\chi$')
set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
ylabel('$\sigma_w / r_\chi$')
set(gca,'YTick', sws,'YTickLabel', ylabels)
fcut = gcf;

% setfig 'b3'
% convergence(convergence > 10) = 10;
% convergence(convergence < 1e-10) = 1e-10;
% contour(B, W, log10(convergence), 'Fill', 'on')
% title('bicgstab convergence (log(estimated error))')
% c = colorbar;
% xlabel(c, 'n = m = 9, N = 3751, \Delta x = \Delta y = 3.33e-4')
% xlabel('\sigma_B / r_\chi')
% set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
% ylabel('\sigma_W / r_\chi')
% set(gca,'YTick', sws,'YTickLabel', ylabels)
% fcon = gcf;
% 
% setfig 'b4'
% contour(B, W, iters, 'Fill', 'on')
% title('number of iterations (max. 100)')
% c = colorbar;
% xlabel(c, 'n = m = 9, N = 3751, \Delta x = \Delta y = 3.33e-4')
% xlabel('\sigma_B / r_\chi')
% set(gca,'XTick', sbs,'XTickLabel', xlabels, 'XTickLabelRotation', 90)
% ylabel('\sigma_W / r_\chi')
% set(gca,'YTick', sws,'YTickLabel', ylabels)
% fiter = gcf;

exportfig(fmid, '../../../images/sigma_depedance_error_gau', '-pdf');
% exportfig(fcut, 'images/sigma_depedance_cuttof_gau', '-pdf');