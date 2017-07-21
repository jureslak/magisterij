prepare

load([plotdatapath 'hertzian_convergence.mat']);

[typenum, simnum, ~] = size(data);

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
Ns = data(1, :, 1);
best = polyfit(log(Ns), log(data(2, :, 2)), 1);
for i = 1:typenum
    plot(Ns, data(i, :, 2), [markers{i}, '-'], 'Color', colors{i})
end
plot(Ns, exp(best(2))*Ns.^best(1), '--k');
text(0.55, 0.32, sprintf('$k = %.2f$', best(1)), 'Units', 'normalized')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$N$')
ylabel('$L_\infty$ napaka')
xlim([1e2, 2e6])
ylim([0.05, 0.5])
set(gca, 'YTick', [0.05, 0.1, 0.2, 0.5])
legend(legendvals)

% f2 = setfig('b2');
% h = area(Ns, time/60);
% for i = 1:6, h(i).FaceColor = colors{i+2}; end
% legend('grajenje domene','ra\v{c}unanje funkcij oblike',...
%        'grajenje matrike', 'ILUT', 're\v{s}evanje sistema',...
%        'ra\v{c}unanje napetosti', 'Location', 'NW');
% %set(gca, 'Yscale', 'log');
% %set(gca, 'Xscale', 'log');
% % set(gca, 'TickDir','out')
% set(gca, 'Layer', 'top')
% xlabel('$N$')
% ylabel('\v{c}as [min]')
% xlim([-inf, inf])


% exportfig(f1, [imagepath 'hertzian_convergence'], '-pdf')
% exportfig(f2, '../../../images/hertzian_time', '-pdf')