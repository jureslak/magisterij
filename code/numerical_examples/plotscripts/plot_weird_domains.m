prepare

load([plotdatapath 'weird_domains.mat'])

f1 = setfig('b1', 'Position', [2000 200 1500 500]);
set(gcf, 'Color', 'w')
s = subplot(1, 3, 2);
daspect([1 1 1])
hold on
axis off
grid off
box off
xint = lpos(2, ltype > 0);
yint = lpos(1, ltype > 0);
xbnd = lpos(2, ltype < 0);
ybnd = lpos(1, ltype < 0);
plot(xint, yint, 'ok')
plot(xbnd, ybnd, 'ok', 'MarkerFaceColor', 'k')
posmid = get(gca, 'Position');
posmid(1) = 0.44;
set(gca, 'Position', posmid);

s = subplot(1, 3, 1);
hold on
axis off
grid off
box off
xint = dpos(1, dtype > 0);
yint = dpos(2, dtype > 0);
xbnd = dpos(1, dtype < 0);
ybnd = dpos(2, dtype < 0);
plot(xint, yint, 'ok')
plot(xbnd, ybnd, 'ok', 'MarkerFaceColor', 'k')
ylim([-inf, inf])
xlim([-inf, inf])
daspect([1 1 1])
pos = get(gca, 'Position');
pos(1) = 0.01;
pos(3) = 0.5;
pos(4) = posmid(4);
set(gca, 'Position', pos);

subplot(1, 3, 3)
daspect([1 1 1])
hold on
axis off
grid off
box off
xint = cpos(1, ctype > 0);
yint = cpos(2, ctype > 0);
zint = cpos(3, ctype > 0);
xbnd = cpos(1, ctype < 0);
ybnd = cpos(2, ctype < 0);
zbnd = cpos(3, ctype < 0);
gray = [0 0 0];
plot3(xint, yint, zint, 'ok')
scatter3(xbnd, ybnd, zbnd, 30, 'o', 'filled',...
    'MarkerFaceColor', gray,...
    'MarkerFaceAlpha', 1)
az = -67.1000;
el = 22.8000;
view(az, el)
pos = get(gca, 'Position');
pos(1) = 0.63;
pos(4) = posmid(4);
pos(3) = 0.3;
set(gca, 'Position', pos);

figure(f1);
exportfig(f1, [imagepath 'weird_domains'], '-pdf')