P = 1000;
L = 10;
D = 5;
E = 72.6e9;
nu = 0.33;

x = linspace(0, L, 100);
y = linspace(-D/2, D/2, 10);
[X, Y] = meshgrid(x, y);
[sxx, syy, sxy, u, v] = cantilever_beam_analytical(X, Y, P, L, D, E, nu);

%%
close all

f1 = setfig('b1');
f = 1e6;
contourf(X+f*u, Y+f*v, sxx, 20);
daspect([1 1 1])
colorbar