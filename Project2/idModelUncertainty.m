function [ymin, ymax] = idModelUncertainty(m, x, horizon, n, alpha)
%Model estimation uncertainty
% input: m, regression model
%        x, regressors
%        horizon, the horizon of history to use
%        n, number of random pertubations
%        alpha, confidence level
% output: ymin, minimum predicted values
%         ymax, maximum predicted values

d = size(m.theta, 1);
thetaRand = randn(d, n);
ys = [];
m_tmp = m; % a temporary regression model (clone of input model)
for i = 1:n 
  thetaRand_i = thetaRand(:, i);
  thetaDelta_i = thetaRand_i/norm(thetaRand_i)*chi2inv(alpha, d);
  m_tmp.theta = m.theta + m.variance^0.5*thetaDelta_i;
  ys = [ys, idpredict(m_tmp, x, horizon)];  % append pred. values with altered theta
end
ymin = min(ys, [], 2); 
ymax = max(ys, [], 2);


end

