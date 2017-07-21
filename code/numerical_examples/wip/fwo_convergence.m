prepare
datafile = [datapath 'fwo_convergence_wip.h5'];
info = h5info(datafile);


a = h5readatt(datafile, '/', 'a');
c = h5readatt(datafile, '/', 'c');
p0 = h5readatt(datafile, '/', 'p0');
R = h5readatt(datafile, '/', 'radius');
COF = h5readatt(datafile, '/', 'COF');
Q = h5readatt(datafile, '/', 'Q');
F = h5readatt(datafile, '/', 'force');
sigmaAxial = h5readatt(datafile, '/', 'sigmaAxial');

anal_approx = 2*p0*sqrt(COF*Q/F) + sigmaAxial;

femnizname = [plotdatapath sprintf('fem_niz/surface_sxx_r%d_mu%g.txt', 1000*R, COF)];
femnizdata = dlmread(femnizname, '\t');

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

f2 = setfig('b2');
% f3 = setfig('b3');
colors = colormap;
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');
    pos = h5read(datafile, [name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);
        vm = von_mises(sxx, syy, sxy);

        displ = h5read(datafile, [name '/displacement']);
        u = displ(1, :);
        v = displ(2, :);

        data(j, i, 1) = N;
        data(j, i, 2) = abs(max(max(sxx)) - anal_approx)/anal_approx;
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');
        
        time(i, :) = [h5readatt(datafile, name, 'time_domain');
                      h5readatt(datafile, name, 'time_shapes');
                      h5readatt(datafile, name, 'time_construct');
                      h5readatt(datafile, name, 'time_lut');
                      h5readatt(datafile, name, 'time_solve');
                      h5readatt(datafile, name, 'time_post')];
    end
%     data(i, 3) = time;

    top = find(y == max(y));
    [~, I] = sort(x(top));
    top = top(I);
    
%     figure(f2)
    plot(x(top)/a, sxx(top), 'o-', 'Color', colors(i, :));
%     
%     [mid, xmid] = closest(pos, 1, 0);
%     figure(f3)
%     plot(sxx(mid)/p0, y(mid)/b, 'o', 'Color', colors(i, :));


    fprintf('point %d/%d \r', i, simnum);
end
plot(femnizdata(:, 1)/a/10^3, femnizdata(:, 2)*10^6, '.', 'Color', 'r')
xlim([-3 3])

%%

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
legendvals = {'G9', 'MON9'};
colors = {
    [1,123,206]/256,
    [255,169,0]/256,
%     [1,56,147]/256,
%     [0,98,199]/256,
    [253,94,91]/256,
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


f1 = setfig('b1');
% best = polyfit(log(Ns), log(data(2, :, 2)), 1);
for i = 1:typenum
    plot(data(1, :, 1), data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
% plot(Ns, exp(best(2))*Ns.^best(1), '--k');
% text(0.55, 0.32, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log')
xlabel('$N$')
ylabel('max $\sigma_{xx}$')
% xlim([1e2, 2e6])
% ylim([0.05, 0.5])
% set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(legendvals)


% f3 = setfig('bo11');
% scatter(x, u);
% scatter(x, v);

% f4 = setfig('b3');
% scatter(x/a, y/a, 5, vm, 'filled');
% daspect([1 1 1])

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


% exportfig(f1, '../../../images/hertzian_convergence', '-pdf')
% exportfig(f2, '../../../images/hertzian_time', '-pdf')