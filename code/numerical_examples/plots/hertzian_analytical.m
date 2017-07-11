function [sxx, syy, sxy] = hertzian_analytical(x, z, b, p0)
% Returns the analytical solution for stress for hertzian contact.
% x in (-inf, inf), z in (-inf, 0].

bxz = b^2 - x.^2 + z.^2;
xz = 4*x.^2.*z.^2;
koren = sqrt(bxz.^2 + xz);
m2 = 1/2 * (koren + bxz);
n2 = 1/2 * (koren - bxz);
m = sqrt(m2);
n = sign(x) .* sqrt(n2);
mpn = m2 + n2;
zmn = (z.^2 + n2) ./ mpn;
sxx = -p0 / b * (m .* (1 + zmn) + 2*z);
syy = -p0 / b * m .* (1 - zmn);
sxy = p0 / b * n .* ((m2 - z.^2) ./ mpn);

end