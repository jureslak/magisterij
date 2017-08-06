prepare
load([plotdatapath 'cantilever_beam_time.mat']);

markers = {'+','o','*','x','s','d','^','v','<','>','p','h'};
legendvals = {'G9 -- displacement', 'G9 -- stress', 'MON9 -- displacement',...
    'MON9 -- stress'};
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

Ns = data(1, :, 1);

f2 = setfig('b2');
h = area(Ns, time/60);
for i = 1:6, h(i).FaceColor = colors{i+2}; end
legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
       'grajenje matrike', 'ILUT', 're\v{s}evanje sistema',...
       'ra\v{c}unanje napetosti', 'Location', 'NW');
%set(gca, 'Yscale', 'log');
%set(gca, 'Xscale', 'log');
% set(gca, 'TickDir','out')
set(gca, 'Layer', 'top')
xlabel('$N$')
ylabel('\v{c}as [min]')
xlim([-inf, inf])

