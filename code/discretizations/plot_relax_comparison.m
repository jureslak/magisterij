clear
close all
randpos = hdf5read('circle_data.h5', '/random/positions');
randtype = hdf5read('circle_data.h5', '/random/types');
randint = find(randtype > 0);
randbnd = find(randtype < 0);

relpos = hdf5read('circle_data.h5', '/relax/positions');
reltype = hdf5read('circle_data.h5', '/relax/types');
relint = find(reltype > 0);
relbnd = find(reltype < 0);

f1 = setfig('b1');
set(gcf, 'Color', 'w')
subplot(1, 2, 2)
daspect([1 1 1])
hold on
axis off
grid off
box off
plot(relpos(1, relint), relpos(2, relint), 'ok')
plot(relpos(1, relbnd), relpos(2, relbnd), 'ok', 'MarkerFaceColor', 'k')

subplot(1, 2, 1)
daspect([1 1 1])
hold on
axis off
grid off
box off
plot(randpos(1, randint), randpos(2, randint), 'ok')
plot(randpos(1, randbnd), randpos(2, randbnd), 'ok', 'MarkerFaceColor', 'k')

% figure(f1)
exportfig(f1, '../../images/relax_comparison', '-pdf')