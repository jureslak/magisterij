prepare
load([plotdatapath 'poisson_square_parallel.mat'])
[typenum, simnum, ~] = size(data);

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
legendvals = cell(typenum, 1);
for i = 1:typenum
    legendvals{i} = sprintf('%d', threads(i));
end
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
    plot(data(i, :, 1), data(1, :, 2)./data(i, :, 2), [markers{i}, '-']);
end

set(gca, 'XScale', 'log')
% title('Speedup in parallel execution using Pardiso solver.')
xlabel('$N$')
ylabel('$t_k$ / $t_1$')
h = legend(legendvals, 'Location', 'NE');
title(h, '$k$', 'interpreter', 'latex')
xlim([2e3, 6e6])

f2 = setfig('b2');
grid minor
thidx = 1;
areadata = reshape(time(thidx, :, :), length(time(thidx, :, :)), []);
h = area(data(thidx, :, 1), areadata/60);
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
grid minor
thidx = 5;
areadata = reshape(time(thidx, :, :), length(time(thidx, :, :)), []);
h = area(data(thidx, :, 1), areadata/60);
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

exportfig(f1, [imagepath 'poisson_square_speedup'], '-pdf');
exportfig(f2, [imagepath 'poisson_square_time_distribution_1'], '-pdf');
exportfig(f3, [imagepath 'poisson_square_time_distribution_12'], '-pdf');