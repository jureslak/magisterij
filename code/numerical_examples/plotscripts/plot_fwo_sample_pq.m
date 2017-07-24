nu = 0.33;
E = 72.1e9;
Estar = E/2/(1-nu^2);
F = 543;
Q = 155;
R = 0.05;
t = 0.004;
mu = 0.85;
sax = 1e8;
p0 = sqrt(F/t*Estar/pi/R);
a = 2*sqrt(F*R/t/pi/Estar);
c = a*sqrt(1-Q/mu/F);
e = a*sax/4/mu/p0;
p = @(x) p0 * sqrt(1-(x/a).^2) .* (abs(x) < a);
q = @(x) mu*p(x) .* (c <= abs(x-e)) + ...
         mu*p0*(sqrt(1-(x/a).^2) - c/a*sqrt(1-((x-e)/c).^2)) .* (abs(x-e) <= c);

xx = linspace(-1.5*a, 1.5*a, 1000);
close all
f1 = setfig('b1');
plot(xx/a, p(xx)/p0, '--', 'Color', [1 1 1]*0.5);
plot(xx/a, q(-xx)/p0, '-k')
legend({'$p(x)$', '$q(x)$'}, 'FontSize', 14)
legend boxoff
axis off
box off

exportfig(f1, [imagepath 'fwo_sample_pq_graph'], '-pdf')