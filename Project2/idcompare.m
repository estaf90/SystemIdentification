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
plot(y_sim, 'LineWidth', 1.5);
grid on;

[y_pred_min, y_pred_max] = idModelUncertainty(m, z, horizon, 100, 0.95);

t = 1:size(z, 1);
patch([t, fliplr(t)], [y_pred_max', fliplr(y_pred_min')], 'r', 'FaceAlpha',0.3, 'EdgeColor' , 'none');

mse_pred = mse(y, y_pred); mse_sim = mse(y, y_sim);
legend('true', ...
sprintf('%d-step prediction (mse=%.1e)', horizon, mse_pred), ...
sprintf('simulated (mse=%.1e)', mse_sim))
title(m.label)

end

