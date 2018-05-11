close all; clear all;

fname = fullfile('data','robot_arm.dat');
data = load(fname);
u = data(:,1); y = data(:,2); t = (1:length(u))';
ut = u(1:end/2); yt = y(1:end/2); tt = t(1:end/2);
uv = u(end/2+1:end); yv = y(end/2+1:end); tv = t(end/2+1:end);
clear fname;

% % Understanding data

figure
subplot(2,1,1); plot(tt, ut, 'b'); hold on; plot(tv,uv,'r'); 
ylabel('torque'); title('input'); legend('train','valid');
subplot(2,1,2); plot(tt, yt, 'b'); hold on; plot(tv, yv,'r'); 
ylabel('acceleration');title('output'); legend('train','valid');

p_utut = periodogram(ut);
p_uvuv = periodogram(uv);
figure
periodogram([ut,uv]); legend('train','valid');

Ts = 1;  % I dont know about this value!
zut = iddata(yt, ut, Ts);
zuv = iddata(yv, uv, Ts);
ge_ut = etfe(zut);
ge_uv = etfe(zuv);
% gs = spa(z1);  % Smoothed spectral estimate
figure; 
bode(ge_ut, ge_uv); legend('train','valid'); grid on;


% auto- and cross-correlation

figure
cv = 2/sqrt(length(u)); %critival value, using normal assumption
subplot(2,2,1); autocorr(u, 'NumLags', 50); title('autocorr(u(t))')
subplot(2,2,3); [r_yu, lags] = xcorr(y, u, 50, 'coeff'); stem(lags, r_yu, '.')
hold on; plot(lags, cv*ones(size(lags)), 'k', lags, -cv*ones(size(lags)), 'k');
title('xcorr(y(t), u(t))')
subplot(2,2,4); autocorr(y, 'NumLags', 50); title('autocorr(y(t))')

% 
% % Data pre-processing
% figure
% subplot(1,2,1); plot(u); title('input'); hold on
% subplot(1,2,2); plot(y); title('output');
% % clear that there are no need of trends
% 
% % Try linear model (wrong)
% m_arx = arx(data, [10,10,0]);
