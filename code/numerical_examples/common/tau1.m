function tau = tau1(sxx, syy, sxy)
  tau = sqrt(1/4*(sxx - syy).^2 + sxy.^2);
end