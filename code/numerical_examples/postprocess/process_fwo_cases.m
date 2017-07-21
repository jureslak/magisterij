prepare
datafile = [datapath 'fwo_cases.h5'];
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

save([plotdatapath 'fwo_cases.mat'], 'data');
