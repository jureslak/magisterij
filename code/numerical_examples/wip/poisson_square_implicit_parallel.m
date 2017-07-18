prepare
datafile = [datapath 'poisson_square_parallel_wip_jarjar.h5'];

info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
threads = {'1', '2', '4', '8', '12', '16', '20', '24'};
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');

    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        data(j, i, 1) = N;
%         data(j, i, 2) = h5readatt(datafile, name, 'thread_num');
        data(j, i, 3) = h5readatt(datafile, name, 'timetotal');

    end

    fprintf('point %d/%d \r', i, simnum);
end

f1 = setfig('b1');
for i = 1:typenum
    plot(data(i, :, 1), data(1, :, 3)./data(i, :, 3), 'x-');
end
set(gca, 'XScale', 'log')
title('Speedup in parallel execution using Pardiso solver.')
xlabel('N')
ylabel('$t_n$ / $t_1$')
legend(threads, 'Location', 'SE')