prepare
load([plotdatapath 'hertzian_matrix_example.mat']);

f1 = setfig('b4');
f = @(x) sign(x).* log10(1+abs(x));
cspy(spfun(f, M))
ch = colorbar;
tck = [-15, -10, -5, 0, 5, 10, 15];
lbl = {'$-10^{15}$', '$-10^{10}$', '$-10^5$', '0', '$10^5$', '$10^{10}$', '$10^{15}$'};
set(ch, 'YTick', tck, 'YTickLabel', lbl)
th = title(sprintf('Neni\\v{c}elnih %.2f\\%%', 100*nnz(M)/numel(M)));
set(th, 'Units', 'normalized', 'Position', [0.55, -0.12], ...
  'VerticalAlignment', 'bottom', ...
  'HorizontalAlignment', 'center');

exportfig(f1, [imagepath 'hertzian_matrix_example.pdf'], '-pdf')