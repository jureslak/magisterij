prepare
datafile = [datapath 'cantilever_convergence_jarjar.h5'];
info = h5info(datafile);

P = h5readatt(datafile, '/', 'P');
D = h5readatt(datafile, '/', 'D');
L = h5readatt(datafile, '/', 'L');
I = h5readatt(datafile, '/', 'I');
E = h5readatt(datafile, '/', 'E');
nu = h5readatt(datafile, '/', 'v');

sb = h5readatt(datafile, '/', 'sigmaB')

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(simnum, 6);

% f2 = setfig('b1');
% f3 = setfig('b3');
% colors = colormap;
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');
    pos = h5read(datafile, [name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    
    [asxx, asyy, asxy, au, av] = cantilever_beam_analytical(x, y, P, L, D, E, nu);
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        syy = stress(2, :);
        sxy = stress(3, :);
        
        displ  = h5read(datafile, [name '/disp']);
        u = displ(1, :);
        v = displ(2, :);
        
        erru = max([max(abs(u - au)), max(abs(v - av))]);
        Mu = max([max(abs(au)), max(abs(av))]);

        errxx = max(max(abs(sxx - asxx)));
        erryy = max(max(abs(syy - asyy)));
        errxy = max(max(abs(sxy - asxy)));
        Ms = max([max(abs(sxx)) max(abs(syy)) max(abs(sxy))]);

        data(j, i, 1) = N;
        data(j, i, 2) = erru / Mu;
        data(j, i, 4) = max([errxx, erryy, errxy]) / Ms;
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


    fprintf('point %d/%d %s\r', i, simnum, name);
end

%%

close all
markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
legendvals = {'G9 -- displacement', 'G9 -- stress', 'MON9 -- displacement',...
    'MON9 -- stress'};
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
Ns = data(1, :, 1);
best = polyfit(log(Ns), log(data(2, :, 2)), 1);
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{2*i}, '-'], 'Color', colors{i})
    plot(Ns, data(i, :, 4), [markers{2*i+1}, '--'], 'Color', colors{i})
end
% plot(Ns, exp(best(2))*Ns.^best(1), '--k');
text(0.50, 0.32, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
% xlim([1e2, 2e6])
% ylim([0.05, 0.5])
% set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(legendvals)

f2 = setfig('b2');
h = area(Ns, time/60);
for i = 1:6, h(i).FaceColor = colors{i+2}; end
legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
       'grajenje matrike', 'ILUT', 're\v{s}evanje sistema',...
       'ra\v{c}unanje napetosti', 'Location', 'NW');
%set(gca, 'Yscale', 'log');
%set(gca, 'Xscale', 'log');
% set(gca, 'TickDir','out')
set(gca, 'Layer', 'top')
xlabel('$N$')
ylabel('\v{c}as [min]')
xlim([-inf, inf])


% exportfig(f1, '../../../images/hertzian_convergence', '-pdf')
% exportfig(f2, '../../../images/hertzian_time', '-pdf')