close all; clear; clc;
dbstop if error;

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

disp('Plotting the Empirical Transfer Function Estimate of training and validation data...')
Ts = 1;  % I dont know about this value!
zt = iddata(yt, ut, Ts);  % training data
zv = iddata(yv, uv, Ts);  % validation data
ge_t = etfe(zt);
ge_v = etfe(zv);
figure; 
bode(ge_t, ge_v); legend('train','valid'); grid on;
set(findall(gcf,'type','line'),'linewidth',2)
pause

% auto- and cross-correlation
%correlations(u, y)

%% Model fitting
model_handles = {@arx, @oe, @n4sid, @ssest, @tfest};
k = 5;
[models, mses, freeParams, akaikePre] = modelOrderVsFreeParams(model_handles, zt, zv, 10, k);

figure;
NameArray = {'LineStyle'};
ValueArray = {'-','-','--',':','-.'}';
subplot(1,3,1); p1 = semilogy(mses,'Linewidth',2); ylabel('mse'); xlabel('model order'); grid on
set(p1,NameArray,ValueArray)
subplot(1,3,2); p2 = plot(freeParams,'Linewidth',2); ylabel('free parameters'); xlabel('model order'); grid on
set(p2,NameArray,ValueArray)
subplot(1,3,3); p3 = plot(akaikePre,'Linewidth',2); ylabel('AkaikePre final prediction error'); xlabel('model order'); grid on
set(p3,NameArray,ValueArray)
legend('arx', 'oe', 'ss (n4sid', 'ss (PEM)', 'tf')
pause

% In discrete-time, a transfer function model has very similar form as the
% OE model, one possible exception is that sometimes by adding feedthrough,
% tf model has parameter b0 in the nominator.


%%
disp('pzmaps...')
pzplots(models, 8)
pause

pzplots(models, 6)
pause

%%
disp('Nyquist plot...')
nyquistplots(models, 8)
pause

nyquistplots(models, 6)
pause

%%
disp('Bode plots...')
bode_plots(models, [6, 8], ge_t)
pause

%%
disp('Simulations')
model_order = 6;
peplot(models, model_order, zv, inf)
pause

disp('Predicions')
model_order = 6;
k = 5;
models_idx = cell2mat(cellfun(@(x) ~isequal(x,@oe), model_handles, 'UniformOutput', false));
peplot(models(:,models_idx), model_order, zv, k)
pause

%%
disp('Residual Analysis')
model_order = 6;
residualAnalysis_plot(models, model_order, zt, 1) % Estimation data?
pause

%%
% ARX
disp('Fitting ARX models')
Opt = arxOptions;
k = 10;
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

disp('Bode plots of the models')
figure
bode(arxL, arxM, arxH)
legend('lower order','medium order', 'high order');

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


disp('Comparing different model structures')
arx_opt = arxH; ss_opt = ss2;
yp_arx = predict(arx_opt, zv, 1);
yp_ss = predict(ss_opt, zv, 1);
mse_arx = immse(yv, yp_arx.y); mse_ss = immse(yv, yp_ss.y);
figure
plot(yv,'LineWidth',1.5); hold on; plot(yp_arx.y); plot(yp_ss.y);
title('Trained models')
legend('y', sprintf('ARX (mse=%.3e)', mse_arx), ...
    sprintf('SS (mse=%.3e)', mse_ss));
pause

figure; pzmap(arx_opt);
figure; pzmap(ss_opt);

%% Functions
% auto- and cross-correlation

function correlations(u_in, y_in)
    figure
    cv = 2/sqrt(length(u_in)); %critival value, using normal assumption
    subplot(2,2,1); autocorr(u_in, 'NumLags', 50); title('autocorr(u(t))')
    subplot(2,2,3); [r_yu, lags] = xcorr(y_in, u_in, 50, 'coeff'); stem(lags, r_yu, '.')
    hold on; plot(lags, cv*ones(size(lags)), 'k', lags, -cv*ones(size(lags)), 'k');
    title('xcorr(y(t), u(t))')
    subplot(2,2,4); autocorr(y_in, 'NumLags', 50); title('autocorr(y(t))')
end


function [trained_models, mses, freeParams, akaikePre] = modelOrderVsFreeParams(models, zt, zv, n, k)
    num_models = length(models);
    mses = zeros(n, num_models);
    akaikePre = zeros(n, num_models);
    freeParams = zeros(n, num_models);
    trained_models = cell(n, num_models);
    for ni=1:n
        for mi = 1:num_models
            model = models{mi};
            mfun = functions(model);
            fprintf('Model order %d for %s\n', ni, mfun.function);
            switch mfun.function
                case 'arx'
                    Opt = arxOptions;
                    ns = [ni ni 1];
                    opt_fun = @() model(zt, ns, Opt);
                case 'oe'
                    Opt = oeOptions;
                    ns = [ni ni 1];
                    Opt.Focus = 'simulation';
                    opt_fun = @() model(zt, ns, Opt);
                case 'n4sid'
                    Opt = n4sidOptions;       
%                     Opt.Display = 'on';
                    ns = ni;
                    opt_fun = @() model(zt, ns, 'Form', 'canonical', Opt);
                case 'ssest'
                    Opt = ssestOptions;       
%                     Opt.Display = 'on';     
                    ns = ni;
                    opt_fun = @() model(zt, ns, 'Form', 'canonical', Opt);
                case 'tfest'
                    Opt = tfestOptions;                   
%                     Opt.Display = 'on';                   
                    Opt.WeightingFilter = [];
                    ns = ni;
                    opt_fun = @() model(zt, ns, Opt, 'Ts', 1);
                otherwise
            end
            
            m = opt_fun();
            trained_models{ni, mi} = m;
            freeParams(ni, mi) = sum(m.Report.Parameters.Free);
            yp = predict(m, zv, k); % k-step prediction
            mse = immse(zv.y, yp.y); 
            mses(ni, mi) = mse;
            akaikePre(ni,mi) = m.Report.Fit.FPE;
        end
    end
end

function bode_plots(models, idxs, ge)
    [nrows, ncols] = size(models);
    rows = 1:nrows;
    for i = 1:ncols
        figure;
        leg = cell(1, length(idxs)+1);
        bode(ge); leg{1,1} = 'train'; hold on;
        for k = 1:length(idxs)
            idx = idxs(k);
            j = rows(idx);
            h = bodeplot(models{j,i});
            showConfidence(h, 0.95);
            leg{1, k+1} = strcat(int2str(idx),' (', models{1,i}.Report.method, ')');
        end
        legend(leg)
        set(findall(gcf,'type','line'),'linewidth',2)
    end
end

function pzplots(models, model_order)
    [nrows, ncols] = size(models);
    rows = 1:nrows;
    M = ceil(sqrt(ncols));
    figure;
    for i = 1:ncols
        subplot(M, ceil(ncols/M), i);
        j = rows(model_order);
        h = iopzplot(models{j,i});
        showConfidence(h, 0.95)
        legend(strcat(int2str(model_order),' (', models{1,i}.Report.method, ')'))
    end
end

function nyquistplots(models, model_order)
    [nrows, ncols] = size(models);
    rows = 1:nrows;
    M = ceil(sqrt(ncols));
    figure;
    for i = 1:ncols
        subplot(M, ceil(ncols/M), i);
        j = rows(model_order);
        h = nyquistplot(models{j,i});
        showConfidence(h, 0.95)
        legend(strcat(int2str(model_order),' (', models{1,i}.Report.method, ')'))
    end
end

function peplot(models, model_order, zv, k)
    [nrows, ncols] = size(models);
    rows = 1:nrows;
    figure;
    leg1 = cell(1, ncols+1);
    leg2 = cell(1, ncols);
    P = zeros(ncols, length(zv.y));
    for i = 1:ncols
        j = rows(model_order);
        pred = predict(models{j,i}, zv, k);
        P(i, :) = pred.y;
        e = zv.y - pred.y;
        subplot(2,1,1); plot(pred.y); hold on;
        subplot(2,1,2); plot(e); hold on;
        leg1{1, i} = models{1,i}.Report.method;
        leg2{1, i} = models{1,i}.Report.method;
    end
    subplot(2,1,1); plot(zv.y, 'k'); hold on; leg1{1,end} = 'valid';
    subplot(2,1,1); legend(leg1);
    subplot(2,1,2); legend(leg2);
    figure
    for i = 1:ncols
        for j = i:ncols
            subplot(ncols, ncols, (i-1)*ncols+j)
            if i == j
                plot(P(i,:) - zv.y')
                title(sprintf('valid - %s',models{1,i}.Report.method))
            else
                plot(P(i,:)-P(j,:))
                title(sprintf('%s - %s',models{1,i}.Report.method,models{1,j}.Report.method))
            end
        ylim([-0.1 0.1])    
            
        end
    end
end

function residualAnalysis_plot(models, model_order, zv, k)
    [nrows, ncols] = size(models);
    rows = 1:nrows;
    figure;
    P = zeros(ncols, length(zv.y));
    for i = 1:ncols
        j = rows(model_order);
        pred = predict(models{j,i}, zv, k);
        P(i, :) = pred.y;
        e = zv.y - pred.y;
        acv = autocorr(e); acv = acv/acv(1)*sqrt(length(e));
        subplot(ncols,1,i); plot(acv,'Linewidth',2); title(strcat('Auto-correlation test ',models{1,i}.Report.method));
    end
    figure;
    for i = 1:ncols
        j = rows(model_order);
        pred = predict(models{j,i}, zv, k);
        P(i, :) = pred.y;
        e = zv.y - pred.y;
        acve = autocorr(e);
        acvu = autocorr(zv.u);
        ccv = crosscorr(zv.u,e); ccv = ccv*sqrt(length(e))/sqrt(acve(1)*acvu(1));
        subplot(ncols,1,i); plot(ccv,'Linewidth',2); title(strcat('Cross-correlation test ',models{1,i}.Report.method));
    end
end

