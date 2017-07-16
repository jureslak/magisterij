prepare
load([plotdatapath 'hertzian_domain_too_small.mat']);

[typenum, simnum, ~] = size(data);

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
colors = {
    [1,123,206]/256, % blue
    [255,169,0]/256, % yellow
    [253,94,91]/256, % red
    [57,172,55]/256, % green
    [115,0,171]/256, % purple
    [100,100,100]/256, % gray
    [115,64,33]/256, % brown
    [223,129,9]/256,
};
names = cell(typenum, 1);
for i = 1:typenum
    names{i} = sprintf('%d\\,b', heights(i)/b);
end

% figure(f2)
% xanal = linspace(-3, 3, 200);
% [asxx, asyy, asxy] = hertzian_analytical(xanal*b, 0*xanal, b, p0);
% avm = von_mises(asxx, asyy, asxy);
% plot(xanal, asxx/p0, '-k')
% xlim([-3, 3])
% figure(f3)
% yanal = linspace(-3, 0, 200);
% [asxx, asyy, asxy] = hertzian_analytical(yanal*0, yanal*b, b, p0);
% plot(asxx/p0, yanal, '-k')
% ylim([-3, 0])

f1 = setfig('b1');
for i = 1:typenum
    plot(data(i, :, 1), data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N \cdot H^2_{max}/H^2$')
ylabel('$L_\infty$ napaka')
xlim([5e3, 3e6])
ylim([0.05, 0.3])
set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
l = legend(names, 'Location', 'SW');
title(l, '$H$')

% f2 = setfig('b2');
% top = find(y == max(y));
% [~, I] = sort(x(top));
% top = top(I);
% 
% plot(x(top)/b, sxx(top)/p0, 'o');
% plot(x(top)/b, asxx(top)/p0, '-');
% 
% plot(x(top)/b, syy(top)/p0, 'o');
% plot(x(top)/b, asyy(top)/p0, '-');
% 
% plot(x(top)/b, sxy(top)/p0, 'o');
% plot(x(top)/b, asxy(top)/p0, '-');
% 
% plot(x(top)/b, sv(top)/p0, 'o');
% plot(x(top)/b, asv(top)/p0, '-');
% 
% legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
%        'von mises analytical', 'Location', 'SE')
% xlim([-3, 3])

% 
% f2 = setfig('b2');
% h = area(Ns, time/60);
% for i = 1:6, h(i).FaceColor = colors{i+2}; end
% legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
%        'grajenje matrike', 'ILUT', 're\v{s}evanje sistema',...
%        'ra\v{c}unanje napetosti', 'Location', 'NW');
% %set(gca, 'Yscale', 'log');
% %set(gca, 'Xscale', 'log');
% % set(gca, 'TickDir','out')
% set(gca, 'Layer', 'top')
% xlabel('$N$')
% ylabel('\v{c}as [min]')
% xlim([-inf, inf])