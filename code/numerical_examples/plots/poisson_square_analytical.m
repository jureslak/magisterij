function z = poisson_square_analytical(x, y)
    z = zeros(size(x));
    for k = 1:2:200
        alfa = 0.5*k*pi;
        z = z - k.^(-3).*pi.^(-3).*sin(k.*pi.*x)...
            .*0.5.*(1-(exp(-2*alfa*y) + exp(-2*alfa*(1-y)))/(1 + exp(-2*alfa)));
            %.*sech(alfa)...
            %.*sinh(alfa*(1-y))...
            %.*sinh(alfa*y);
    end
    z = 8*z;
end

