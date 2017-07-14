function [sxx, syy, sxy, u, v] = cantilever_beam_analytical(x, y, P, L, D, E, nu)
% returns analytical solution to timoshenko cantilever beam
% P = loading force
% L = length
% D = height
% E = Young's modulus
% nu = Poisson ratio
I = D^3 / 12;
sxx = I.^(-1).*P.*x.*y;
syy = zeros(size(x));
sxy = (1/2).*I.^(-1).*P.*((1/4).*D.^2+(-1).*y.^2);

u = (1/24).*E.^(-1).*I.^(-1).*P.*y.*(3.*D.^2.*(1+nu)-...
    4.*(3.*L.^2+(-3).*x.^2+(2+nu).*y.^2));
v = (-1/24).*E.^(-1).*I.^(-1).*P.*(3.*D.^2.*(1+nu).*(L+(-1).*x)+...
    4.*(2.*L.^3+(-3).*L.^2.*x+x.^3+3.*nu.*x.*y.^2));
end