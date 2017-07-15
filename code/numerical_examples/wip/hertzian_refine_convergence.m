prepare
datafile = [datapath 'hertzian_refined_convergence_jarjar.h5'];
info = h5info(datafile);

b = h5readatt(datafile, '/', 'b');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

% f2 = setfig('b1');
% f3 = setfig('b3');
% colors = colormap;
for i = 1:simnum
    if i > 22, break, end
    if i < 9, continue, end
   
    for j = 1:typenum
        if i == 22 && j > 5, continue, end
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        N = h5readatt(datafile, name, 'N');
        pos = h5read(datafile, [name '/pos']);
        x = pos(1, :);
        y = pos(2, :);

        testidx = find(x.^2 + y.^2 < (1000*b)^2);% & y.^2 > (b/2)^2);
        [asxx, asyy, asxy] = hertzian_analytical(x, y, b, p0);
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);

        errxx = max(max(abs(sxx(testidx) - asxx(testidx))));
        erryy = max(max(abs(syy(testidx) - asyy(testidx))));
        errxy = max(max(abs(sxy(testidx) - asxy(testidx))));

        data(j, i, 1) = N;
        data(j, i, 2) = max([errxx, erryy, errxy])/p0;
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');
        
        time(i, :) = [h5readatt(datafile, name, 'time_domain');
                      h5readatt(datafile, name, 'time_shapes');
                      h5readatt(datafile, name, 'time_construct');
                      h5readatt(datafile, name, 'time_lut');
                      h5readatt(datafile, name, 'time_solve');
                      h5readatt(datafile, name, 'time_post')];
    end
%     data(i, 3) = time;

%     top = y == max(y);
%     figure(f2)
%     plot(x(top)/b, sxx(top)/p0, 'o', 'Color', colors(i, :));
%     
%     [mid, xmid] = closest(pos, 1, 0);
%     figure(f3)
%     plot(sxx(mid)/p0, y(mid)/b, 'o', 'Color', colors(i, :));


    fprintf('point %d/%d \r', i, simnum);
end

%%

close all
markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
lines = {'-', '--', '-', '--', '-', '--', '-'};
legendvals = {'refine0', 'refine1', 'refine2', 'refine3', 'refine4', 'refine5', 'refine6'};
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
% Ns = data(1, :, 1);
% best = polyfit(log(Ns), log(data(1, :, 2)), 1);
for i = 1:typenum
    plot(data(i, :, 1), data(i, :, 2), [markers{i}, lines{i}], 'Color', colors{i})
end
% plot(Ns, exp(best(2))*Ns.^best(1), '--k');
% text(0.55, 0.32, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([8e3, 1.5e6])
ylim([-inf, 6e-2])
set(gca, 'YTick', [5e-3, 1e-2, 5e-2])
legend(legendvals)

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

% exportfig(f1, '../../../images/hertzian_refined_convergence', '-pdf')
% exportfig(f2, '../../../images/hertzian_refined_time', '-pdf')