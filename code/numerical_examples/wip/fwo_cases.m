prepare
datafile = [datapath 'fwo_cases_wip.h5'];
info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = cell(typenum, simnum, 3);

for i = 1:simnum
    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
        
        N = h5readatt(datafile, name, 'N');
        pos = h5read(datafile, [name '/pos']);
        x = pos(1, :);
        y = pos(2, :);
        
        top = find(y == max(y));
        [~, I] = sort(x(top));
        top = top(I);
    
        stress = h5read(datafile, [name '/stress']);
        sxx = stress(1, :);
        
        a = h5readatt(datafile, name, 'a');
        R = h5readatt(datafile, name, 'R');
        COF = h5readatt(datafile, name, 'COF');
        
        femnizname = [plotdatapath sprintf('fem_niz/surface_sxx_r%d_mu%g.txt', 1000*R, COF)];
        femnizdata = dlmread(femnizname, '\t');

        data{j}{i}{1} = x(top);
        data{j}{i}{2} = sxx(top) / 10^6;
        data{j}{i}{3} = femnizdata;
        data{j}{i}{4} = [a, COF, R];
    end
    fprintf('entry %d/%d \r', i, simnum);
end

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