prepare
datafile = [datapath 'poisson_square_implicit.h5'];
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
        % if strcmp(name, '/mon9/calc0408'), break, end
    end
    fprintf('point %d/%d \r', i, simnum);
end

names = cell(typenum, 1);
for i = 1:typenum
    names{i} = info.Groups(i).Name;
end
save([plotdatapath 'poisson_square_implicit.mat'], 'data', 'names', 'typenum', 'simnum')