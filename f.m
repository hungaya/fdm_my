function f= f(u)
  EPSILON = 1e-6;
  u = abs(u);
  if u < EPSILON
    f = 0;
  else
    f = u*log(u);
  end
end
