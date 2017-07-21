prepare
load([plotdatapath 'fwo_cases.mat']);

[typenum, simnum, ~] = size(data);

f1 = setfig('b1');
for i = 1:typenum
    for j = 1:simnum
        subplot(2, 2, 2*(i-1)+j);
        hold on
        grid on
        box on
        a = data{i}{j}{4}(1);
        COF = data{i}{j}{4}(2);
        R = data{i}{j}{4}(3);
        plot(data{i}{j}{1}/a, data{i}{j}{2}, '.-');    
        plot(data{i}{j}{3}(:, 1)/a/10^3, data{i}{j}{3}(:, 2), '.-');
        title(sprintf('$\\mu$ = %g, $R$ = %g\\,mm', COF, 1000*R));
        legend('MLSM', 'FEM', 'Location', 'NW')
        xlabel('$x/a$')
        ylabel('$\sigma_{xx}$ [MPa]')
        xlim([-2, 2])
    end
end

exportfig(f1, [imagepath 'fwo_cases'], '-pdf');