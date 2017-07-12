format compact
clear
datafile = '../data/convergence_test_1D_rbf1.h5';
info = h5info(datafile);

correct = @(x) cos(1)*x - sin(x);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);

for i = 1:simnum
    pos = h5read(datafile, [info.Groups(1).Groups(i).Name '/pos']);
    anal = correct(pos);
    
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        pos = h5read(datafile, [name '/pos']);
        anal = correct(pos);
        if (isrow(anal)), anal = anal'; end

        sol = h5read(datafile, [name '/sol']);
        N = h5readatt(datafile, name, 'N');
        time = h5readatt(datafile, name, 'timetotal');
       
        M = spconvert(h5read(datafile, [name '/M'])');
        rhs = h5read(datafile, [name '/rhs']);
        % sol = M \ rhs;

        err = max(max(abs(sol - anal)));
        data(j, i, 1) = N;
        data(j, i, 2) = err;
        data(j, i, 3) = h5readatt(datafile, name, 'cutoff');

    %     timem(i) = timefdm;
    %     times(i) = time;
    %     timespart(i) = timepart;
    %     
    %     errs(i) = max(abs(correct(xx)' - yy));
    %     errsfdm(i) = max(abs(correct(linspace(0, 1, length(xx)))' - yyfdm));
    %     Ns(i) = length(xx);

        %if i > 50, break, end
        % if strcmp(name, '/mon9/calc0408'), break, end
        if N == 1024,
            fprintf('\n%s:\n', name);
            shape = M(5, :)
        end
    end
%     fprintf('point %d/%d \r', i, simnum);
end

%%

close all
% keys = {'/gau3','/gau5','/imq3','/imq5','/mon3','/mon5','/mq3','/mq5'};
% vals = {'G3','G5','IMQ3','IMQ5', 'MON3', 'MON5', 'MQ3','MQ5'};
% namemap = containers.Map(keys, vals);
legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = strrep(info.Groups(i).Name, '_', '\_');
end

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
colors = {
    [1,123,206]/256,
    [1,56,147]/256,
    [253,94,91]/256,
    [137,1,1]/256,
    [255,207,0]/256,
    [255,169,0]/256,
    [57,172,55]/256,
    [15,85,48]/256,
};

Ns = data(1, :, 1);

%ok = find(Ns < 1e3);
%best = polyfit(log(Ns), log(data(8, :, 2)), 1);
f1 = setfig('b1');
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
%text(0.4, 0.43, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
% ylim([1e-7, 1e2])
xlim([-inf, 1e5])
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
ylabel('avg. cutoff')