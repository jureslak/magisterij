clear
datafile = '../data/poisson_square_implicit_basis.h5';
info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);

for i = 1:simnum
    pos = h5read(datafile, [info.Groups(1).Groups(i).Name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    anal = poisson_square_analytical(x, y)';
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        sol1 = double(h5read(datafile, [name '/sol']));
        N = h5readatt(datafile, name, 'N');
        time = h5readatt(datafile, name, 'timetotal');

        M = spconvert(h5read(datafile, [name '/M'])');
        rhs = h5read(datafile, [name '/rhs']);
        sol = M \ rhs;
        
        fprintf('Diff at %s = %g\n', name, norm(sol1 - sol));
        
        err = max(max(abs(sol - anal)));
        cutoff = h5read(datafile, [name, '/cutoff']);
        cutoff = mean(reshape(cutoff, [2, length(cutoff)/2]));

        
        assert(norm(M*sol - rhs) < 1e-6,...
            'Error is %f at %s.', norm(M*sol - rhs), name);
        
        data(j, i, 1) = N;
        data(j, i, 2) = err;
        data(j, i, 3) = mean(cutoff);
    end
    fprintf('point %d/%d \r', i, simnum);
end

%%

close all
legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = info.Groups(i).Name;
end

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
colors = {
    [1,56,147]/256,
    [1,123,206]/256,
    [0,98,199]/256,
    [137,1,1]/256,
    [253,94,91]/256,
    [238,16,31]/256,
    [255,207,0]/256,
    [255,169,0]/256,
    [223,129,9]/256,
    [15,85,48]/256,
    [57,172,55]/256,
    [19,131,49]/256
};

Ns = data(1, :, 1);

f1 = setfig('b1');
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([40, 10^5])
legend(legendvals)

f2 = setfig('b2', 'Visible', 'on');
hold off
hold on
for i = 1:typenum
    plot(Ns, data(i, :, 3), [markers{i}, '-'], 'Color', colors{i})
end
set(gca, 'XScale', 'log')%, 'YScale', 'log')
xlim([40, 10^5])
legend(legendvals, 'Location', 'SE')
xlabel('$N$')
ylabel('\v{c}as [$s$]')

f3 = setfig('b3');
name = '/gau_s00300.00/calc0061';
cutoff = h5read(datafile, [name, '/cutoff']);
cutoff = mean(reshape(cutoff, [2, length(cutoff)/2]));

pos = h5read(datafile, [name '/pos']);
sol = h5read(datafile, [name '/sol']);
x = pos(1, :);
y = pos(2, :);
scatter(x, y, 25, sol, 'filled');
colorbar

% exportfig(f1, '../../../images/poisson_square_convergence_basis', '-pdf')
% exportfig(f2, '../../../images/poisson_square_time', '-pdf')