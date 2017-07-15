prepare 

datafile = [datapath 'hertzian_domain_too_small.h5'];
info = h5info(datafile);

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');
max_height = h5readatt(datafile, '/', 'max_height');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

% f2 = setfig('b1');
% f3 = setfig('b3');
% colors = colormap;
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name
        
        N = h5readatt(datafile, name, 'N');
        height = h5readatt(datafile, name, 'height');
        pos = h5read(datafile, [name '/pos']);
        x = pos(1, :);
        y = pos(2, :);

        testidx = find(x.^2 + y.^2 < (10000*b)^2);
        [asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);
        asv = von_mises(asxx, asyy, asxy);
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);
        sv = von_mises(sxx, syy, sxy);

        errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
        erryy = max(max(abs(syy(testidx) - asyy(testidx))));
        errxy = max(max(abs(sxy(testidx) - asxy(testidx))));

        data(j, i, 1) = N * (max_height*b / height)^2;
        data(j, i, 2) = max([errxx, erryy, errxy])/p0;
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');
        
    end
%     data(i, 3) = time;

%     top = y == max(y);
%     figure(f2)
%     plot(x(top)/b, sxx(top)/p0, 'o', 'Color', colors(i, :));
%     
%     [mid, xmid] = closest(pos, 1, 0);
%     figure(f3)
%     plot(sxx(mid)/p0, y(mid)/b, 'o', 'Color', colors(i, :));


    fprintf('point %d/%d %s\r', i, simnum, name);
end

%%
close all
markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
% legendvals = {'MON9'};
colors = {
    [1,123,206]/256,
    [255,169,0]/256,
%     [1,56,147]/256,
%     [0,98,199]/256,
    [253,94,91]/256,
    [19,131,49]/256,
    [255,169,0]/256,
%     [137,1,1]/256,
%     [238,16,31]/256,
    [15,85,48]/256+0.1,
    [57,172,55]/256,
    [1,123,206]/256,
    [255,207,0]/256,
    [223,129,9]/256,
    [19,131,49]/256
};
names = cell(typenum, 1);
for i = 1:typenum
    names{i} = info.Groups(i).Name;
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
xlabel('$N \cdot A/A_{max}$')
ylabel('$L_\infty$ napaka')
% xlim([1e2, 2e6])
% ylim([0.05, 0.5])
set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(names)

f2 = setfig('b2');
top = find(y == max(y));
[~, I] = sort(x(top));
top = top(I);

plot(x(top)/b, sxx(top)/p0, 'o');
plot(x(top)/b, asxx(top)/p0, '-');

plot(x(top)/b, syy(top)/p0, 'o');
plot(x(top)/b, asyy(top)/p0, '-');

plot(x(top)/b, sxy(top)/p0, 'o');
plot(x(top)/b, asxy(top)/p0, '-');

plot(x(top)/b, sv(top)/p0, 'o');
plot(x(top)/b, asv(top)/p0, '-');

legend('sxx', 'asxx', 'syy', 'asyy', 'sxy', 'asxy', 'von mises',...
       'von mises analytical', 'Location', 'SE')
xlim([-3, 3])

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
