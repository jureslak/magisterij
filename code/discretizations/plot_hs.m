close all
Ns = hdf5read('hs_data.h5', '/circle/N');
hs = hdf5read('hs_data.h5', '/circle/h');
Ss = hdf5read('hs_data.h5', '/circle/S');

Ns_relax = hdf5read('hs_data.h5', '/circle_refine/N');
hs_relax = hdf5read('hs_data.h5', '/circle_refine/h');
Ss_relax = hdf5read('hs_data.h5', '/circle_refine/S');


f1 = setfig('b1');
c1 = plot(Ns, hs, 's', Ns_relax, hs_relax, 'o');
c2 = plot(Ns, Ss, '^');
plot(Ns_relax, Ss_relax, 'v');
set(gca, 'YScale', 'log')
xlabel('$N$', 'FontSize', 14)
xlim([0 3000])
leg = legend('$h$', '$h_{\textrm{spro\v{s}\v{c}ena}}$',...
    '$S$', '$S_{\textrm{spro\v{s}\v{c}ena}}$');
leg.FontSize = 14;
uistack(c1, 'top')
uistack(c2, 'top')
set(gca,'Position', [.06 .1 .91 .87])
set(gca,'FontSize', 14)

figure(f1)
exportfig(f1, '../../images/relax_improvement', '-pdf')