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


p_uu = periodogram(u);
p_utut = periodogram(ut);
p_uvuv = periodogram(uv);

Ts = 1;  % I dont know about this value!
zu = iddata(y, u, Ts);
zut = iddata(yt, ut, Ts);
zuv = iddata(yv, uv, Ts);

ge_u = etfe(zu);  % Empirical Transfer Function Estimate
ge_ut = etfe(zut);
ge_uv = etfe(zuv);
% gs = spa(z1);  % Smoothed spectral estimate
figure; 
bode(ge_ut, ge_uv); legend('train','valid'); grid on;
% 
% % Data pre-processing
% figure
% subplot(1,2,1); plot(u); title('input'); hold on
% subplot(1,2,2); plot(y); title('output');
% % clear that there are no need of trends
% 
% % Try linear model (wrong)
% m_arx = arx(data, [10,10,0]);
