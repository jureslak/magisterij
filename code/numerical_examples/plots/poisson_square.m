close all
clear
datafile = '../data/poisson_square.h5';
info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);

for i = 1:simnum
    pos = h5read(datafile, [info.Groups(1).Groups(i).Name '/pos']);
    x = pos(1, :);
    y = pos(2, :);
    anal = poisson_square_analytical(x, y);
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        sol = h5read(datafile, [name '/sol'])';
        N = h5readatt(datafile, name, 'N');
        time = h5readatt(datafile, name, 'timetotal');

        err = max(max(abs(sol - anal)));
        data(j, i, 1) = N;
        data(j, i, 2) = err;
        data(j, i, 3) = time;



    %     timem(i) = timefdm;
    %     times(i) = time;
    %     timespart(i) = timepart;
    %     
    %     errs(i) = max(abs(correct(xx)' - yy));
    %     errsfdm(i) = max(abs(correct(linspace(0, 1, length(xx)))' - yyfdm));
    %     Ns(i) = length(xx);

        %if i > 50, break, end
    end
    fprintf('point %d/%d \r', i, simnum);
end

%%

close all
keys = {'/gau13','/gau5','/gau9','/imq13','/imq5','/imq9','/mon5','/mon9','/mq13','/mq5','/mq9', '/mon6'};
vals = {'G13','G5','G9','IMQ13','IMQ5','IMQ9', 'MON5', 'MON9', 'MQ13','MQ5','MQ9', 'MON6'};
namemap = containers.Map(keys, vals);
legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = namemap(info.Groups(i).Name);
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
xlim([40, Inf])
legend(legendvals)
figure(f1)

f2 = setfig('b2');
for i = 1:typenum
    plot(Ns, data(i, :, 3), [markers{i}, '-'], 'Color', colors{i})
end
% set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([40, Inf])
legend(legendvals, 'Location', 'SE')
xlabel('$N$')
ylabel('\v{c}as [$s$]')
figure(f2)

%exportfig(f1, '../../../images/poisson_square_convergence', '-pdf')
%exportfig(f2, '../../../images/poisson_square_time', '-pdf')