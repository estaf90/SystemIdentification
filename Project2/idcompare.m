function idcompare(m, z, horizon)
%IDCOMPARE plot the predicted/simulated output with true output
%   Detailed explanation goes here

% [n, d] = size(z);
y = z(:,1);

y_pred = idpredict(m, z, horizon);
y_sim = idsimulate(m, z);
figure; hold on;
plot(y, 'LineWidth', 1.5);
plot(y_pred, 'LineWidth', 1.5);
if strcmp(m.model, 'ARX')
    plot(y_sim, 'LineWidth', 1.5);
    
    ylim([min([min(y_pred), min(y_sim), min(y)]), ...
         max([max(y_pred), max(y_sim), max(y)])])
else
    ylim([min([min(y_pred), min(y)]), ...
         max([max(y_pred), max(y)])])
end

grid on;

[y_pred_min, y_pred_max] = idModelUncertainty(m, z, horizon, 100, 0.95);

t = 1:size(z, 1);
patch([t, fliplr(t)], [y_pred_max', fliplr(y_pred_min')], 'r', 'FaceAlpha',0.3, 'EdgeColor' , 'none');

mse_pred = mse(y, y_pred); mse_sim = mse(y, y_sim);
if strcmp(m.model, 'OE')
    h_leg = legend('true', ...
    sprintf('Prediction (mse=%.1e)', mse_pred));
else
    h_leg = legend('true', ...
    sprintf('%d-step prediction (mse=%.1e)', horizon, mse_pred), ...
    sprintf('simulated (mse=%.1e)', mse_sim));
end
set(h_leg, 'Location', 'Best')
title(m.label)

% sys = idpoly([1, m.A'], [zeros(1, m.n(3)), m.B']);
% sys = setcov(sys, m.variance);
% 
% figure
% h = iopzplot(sys);
% showConfidence(h, 2);
% 
% figure
% h = bodeplot(sys);
% showConfidence(h, 2);
% 
% figure()
% h = nyquistplot(sys);
% showConfidence(h, 2);

end

