prepare
datafile = [datapath 'fwo_table_wip.h5'];
info = h5info(datafile);

radnum = length(info.Groups);
cofnum = length(info.Groups(1).Groups);
simnum = length(info.Groups(1).Groups(1).Groups);

cofs = [0.3 0.85 2];
rads = [0.01, 0.05];

data = zeros(radnum, cofnum, simnum, 3);

for i = 1:radnum
    for j = 1:cofnum
        for k = 1:simnum
%             if k > 4, break, end
            grp = info.Groups(i).Groups(j).Groups(k);
            name = grp.Name;
            
            a = h5readatt(datafile, name, 'a');
            pos = h5read(datafile, [name '/pos']);
            x = pos(1, :);
            y = pos(2, :);
            n = sort(sqrt((x+a).^2 + y.^2));
            md = n(2)*10^6
            
            

            N = h5readatt(datafile, name, 'N');
            stress = h5read(datafile, [name '/stress']);
            maxsxx = max(stress(1, :));

            data(i, j, k, 1) = N;
            data(i, j, k, 2) = maxsxx;
            
            fprintf('entry %d/%d \r', k, simnum);

        end
    end
end

%%

% table = data(:, :, end-2, 2)/10^6
% table = data(:, :, end-1, 2)/10^6
table = data(:, :, :, 2)/10^6
