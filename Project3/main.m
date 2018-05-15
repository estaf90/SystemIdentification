close all; clear all;

fname = fullfile('data','robot_arm.dat');
disp('Loading data...')
data = load(fname);
u = data(:,1); y = data(:,2); t = (1:length(u))';
sf = 2;  % split factor
ut = u(1:end*(1-1/sf)); yt = y(1:end*(1-1/sf)); tt = t(1:end*(1-1/sf));
uv = u(end*(1-1/sf)+1:end); yv = y(end*(1-1/sf)+1:end); tv = t(end*(1-1/sf)+1:end);
clear fname;

% % Understanding data
disp('Understanding the data')
disp('Plotting the input and the output signals...')
figure
subplot(2,1,1); plot(tt, ut, 'b'); hold on; plot(tv,uv,'r'); 
ylabel('torque'); title('input'); legend('train','valid');
subplot(2,1,2); plot(tt, yt, 'b'); hold on; plot(tv, yv,'r'); 
ylabel('acceleration');title('output'); legend('train','valid');

pause

% disp('Plotting periodograms of training and validation data...')
% p_utut = periodogram(ut);
% p_uvuv = periodogram(uv);
% figure
% periodogram([ut,uv]); legend('train','valid');
% pause

disp('Plotting the Empirical Transfer Function Estimate of training and validation data...')
Ts = 1;  % I dont know about this value!
zt = iddata(yt, ut, Ts);  % training data
zv = iddata(yv, uv, Ts);  % validation data
ge_t = etfe(zt);
ge_v = etfe(zv);
% gs = spa(z1);  % Smoothed spectral estimate
figure; 
bode(ge_t, ge_v); legend('train','valid'); grid on;
pause

% auto- and cross-correlation

figure
cv = 2/sqrt(length(u)); %critival value, using normal assumption
subplot(2,2,1); autocorr(u, 'NumLags', 50); title('autocorr(u(t))')
subplot(2,2,3); [r_yu, lags] = xcorr(y, u, 50, 'coeff'); stem(lags, r_yu, '.')
hold on; plot(lags, cv*ones(size(lags)), 'k', lags, -cv*ones(size(lags)), 'k');
title('xcorr(y(t), u(t))')
subplot(2,2,4); autocorr(y, 'NumLags', 50); title('autocorr(y(t))')


% Model fitting
disp('Fitting ARX models')
Opt = arxOptions;
k = 5;
arxL = arx(zt, [4 4 1], Opt);
arxM = arx(zt,[6 6 1], Opt);
arxH = arx(zt,[8 8 1], Opt);
predL = predict(arxL, zv, k);
predM = predict(arxM, zv, k);
predH = predict(arxH, zv, k);
mseH = immse(yv, predH.y); mseM = immse(yv, predM.y); mseL = immse(yv, predL.y);
figure
plot(yv,'LineWidth',1.5); hold on; plot(predL.y); plot(predM.y); plot(predH.y);
title('Trained ARX models')
legend('y', sprintf('lower order (%.3e)', mseL),...
    sprintf('medium order (%.3e)', mseM), ...
    sprintf('high order (%.3e)', mseH));
pause

disp('The low order model is to simple')
figure
plot(yv,'LineWidth',1.5); hold on; plot(predM.y); plot(predH.y);
title('Trained ARX models')
legend('y', sprintf('medium order (%.3e)', mseM), ...
    sprintf('high order (%.3e)', mseH));
disp('The training and validation data include very little variation which makes the model estimation prone to overfit without showing in the mse (variance part)!')
disp('We need to be careful that the model is not overfitting')
pause

% State space model estimation 
disp('Estimating State Space models')
                      
Options = n4sidOptions;                             
Options.Display = 'on';
ss1 = n4sid(zv, 8, Options);
ss2 = n4sid(zv, 10, 'Form', 'canonical', Options);
pred_ss1 = predict(ss1, zv, k);
pred_ss2 = predict(ss2, zv, k);
mse_ss1 = immse(yv, pred_ss1.y); mse_ss2 = immse(yv, pred_ss2.y);
figure
plot(yv,'LineWidth',1.5); hold on; plot(pred_ss1.y); plot(pred_ss2.y);
title('Trained State Space models')
legend('y', sprintf('Free (mse=%.3e)', mse_ss1), ...
    sprintf('OC (mse=%.3e)', mse_ss2));

pause