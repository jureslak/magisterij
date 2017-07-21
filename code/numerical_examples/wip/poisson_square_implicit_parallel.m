prepare
datafile = [datapath 'poisson_square_parallel_wip_jarjar.h5'];

info = h5info(datafile);

typenum = length(info.Groups);
simnum = length(info.Groups(1).Groups);

data = zeros(typenum, simnum, 3);
time = zeros(typenum, simnum, 5);
for i = 1:simnum
    name = info.Groups(1).Groups(i).Name;

    N = h5readatt(datafile, name, 'N');

    for j = 1:typenum
        grp = info.Groups(j).Groups(i);
        name = grp.Name;
    
        data(j, i, 1) = N;
%         data(j, i, 2) = h5readatt(datafile, name, 'thread_num');
        data(j, i, 3) = h5readatt(datafile, name, 'time_total');

        
        time(j, i, :) = [h5readatt(datafile, name, 'time_domain');
                         h5readatt(datafile, name, 'time_shapes');
                         h5readatt(datafile, name, 'time_construct')
                         h5readatt(datafile, name, 'time_compute')
                         h5readatt(datafile, name, 'time_solve');];
    end

    fprintf('point %d/%d \r', i, simnum);
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
    plot(data(i, :, 1), data(1, :, 3)./data(i, :, 3), [markers{i}, '-']);
end

set(gca, 'XScale', 'log')
% title('Speedup in parallel execution using Pardiso solver.')
xlabel('N')
ylabel('$t_n$ / $t_1$')
threads = {'1', '2', '4', '8', '12', '16', '20', '24'};
legend(threads, 'Location', 'NE')

f2 = setfig('b2');
areadata = reshape(time(1, :, :), length(time(1, :, :)), []);
h = area(data(1, :, 1), areadata/60);
for i = 1:5, h(i).FaceColor = colors{i+2}; end
legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
       'grajenje matrike', 'ra\v{c}unanje razcepa', 're\v{s}evanje sistema',...
       'Location', 'NW');
%set(gca, 'Yscale', 'log');
%set(gca, 'Xscale', 'log');
% set(gca, 'TickDir','out')
set(gca, 'Layer', 'top')
xlabel('$N$')
ylabel('\v{c}as [min]')
xlim([-inf, inf])

f3 = setfig('b3');
areadata = reshape(time(5, :, :), length(time(5, :, :)), []);
h = area(data(5, :, 1), areadata/60);
for i = 1:5, h(i).FaceColor = colors{i+2}; end
legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
       'grajenje matrike', 'ra\v{c}unanje razcepa', 're\v{s}evanje sistema',...
       'Location', 'NW');
%set(gca, 'Yscale', 'log');
%set(gca, 'Xscale', 'log');
% set(gca, 'TickDir','out')
set(gca, 'Layer', 'top')
xlabel('$N$')
ylabel('\v{c}as [min]')
xlim([-inf, inf])