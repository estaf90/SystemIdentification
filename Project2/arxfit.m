function m = arxfit(Z, N)
  na = N(1); % number of past outputs in regressor
  nb = N(2); % number of past inputs
  nk = N(3); % delays in input signal
  
  [X, Y] = arxinput(Z, [na nb nk]);
  
  m = LinRegress(X, Y);
  m.model = 'ARX';
  m.label = sprintf('ARX [na=%d nb=%d nk=%d]', na,nb,nk);
  m.n = [na, nb, nk];
  m.A = m.theta(1:na);
  m.B = m.theta((na+1):end);
end