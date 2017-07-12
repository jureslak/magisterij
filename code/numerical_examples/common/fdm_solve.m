%function s = fdm_solve(f, a, b, N)

f = @sin;
a = 0;
b = 0;
N = 101;

M = zeros(N, N);
rhs = zeros(N, 1);
h = 1/(N-1);
for i=3:N-2
    M(i, i-2:i+2) = [-1/12, 4/3, -5/2, 4/3, -1/12]/h^2;
    %M(i, i-1:i+1) = [1 -2 1]/h^2;
    rhs(i) = f((i-1)*h);
end
M(2, 1:5) = [-20 11 6 4 -1]/h^2/12;
M(end-1, end-4:end) = [-1 4 6 11 -20]/h^2/12;
rhs(2) = f(h);
rhs(end-1) = f(1-h);

rhs(1) = a;
M(1, 1) = 1;

rhs(end) = b;
% M(end, end) = 1;
M(end, end-4:end) = [1/4, -4/3,  3, -4, 25/12]/h;
% M(end, end-2:end) = [1/2 -2 3/2]/h;

s = M \ rhs;

close all
setfig('b1');
plot(linspace(0, 1, N), s, 'o-')
xa = linspace(0, 1, 10000);
plot(xa, cos(1)*xa-sin(xa), '--')
plot(xa, sin(1)*xa-sin(xa), '--')

% end