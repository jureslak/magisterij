function v = von_mises(sxx, syy, sxy)
  v = sqrt(sxx.^2 - sxx.*syy + syy.^2 + 3*sxy.^2);
end
