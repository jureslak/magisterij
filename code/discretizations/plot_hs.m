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
xlabel('N')
xlim([0 3000])
legend('$h$', '$h_{relax}$', '$S$', '$S_{relax}$')
uistack(c1, 'top')
uistack(c2, 'top')
set(gca,'Position', [.06 .1 .91 .87])

figure(f1)
exportfig(f1, '../../images/relax_improvement', '-pdf')