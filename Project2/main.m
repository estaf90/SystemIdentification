close all; clear all;
dbstop if error
% 1.a
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load control
end

x = gen_input(100, 'white');

a = [0.2, -0.3]; na = length(a);
b = [0, 0.4, -0.2]; nb = length(b)-1;
nk = 1;
nfit = [na, nb, nk];
y = gen_arx_data(x, a, b, 0);
yv = gen_arx_data(x, a, b, 1e-2);

Z = [y, x];

m_arx = arxfit(Z, nfit);

disp('1 (a): ARX expected: [0.2, -0.3, 0.4, -0.2]')
disp(m_arx.theta')  % Verify that the parameters are found
idcompare(m_arx, Z, 1)

y = gen_arx_data(x, a, b, 1e-1);
yvalid = gen_arx_data(x, a, b, 1e-1);
Z = [y, x];
Zvalid = [yvalid, x];
m_arx = arxfit(Z, nfit);

% 1.b
H = id2tf(m_arx);
% ltiview({'step'}, H);

%1.c
disp('1 (c)')

%1.d
disp('1 (d)')
idcompare(m_arx, Zvalid, 1)
idcompare(m_arx, Zvalid, 5)
idcompare(m_arx, Zvalid, 10)
% idcompare(m_arx, Zvalid, 50) % should be the same as simulation

%2
disp('2')
noise = 5e-2;
n = 500;
x = gen_input(n, 'steps');
y_oe = gen_oe_data(x, a, b, noise);
Z_oe = [y_oe, x];

oe_fit = oe(Z_oe, [nb, na, nk]);

m2_oe = oefit(Z_oe, nfit);
m2_arx = arxfit(Z_oe, nfit);

x = gen_input(n, 'steps');
y_oe = gen_oe_data(x, a, b, noise);
Z_oe = [y_oe, x];
horizon = 1;
idcompare(m2_oe, Z_oe, horizon)
idcompare(m2_arx, Z_oe, horizon)


n = 1e4;

x = gen_input(n, 'steps');
y_oe = gen_oe_data(x, a, b, noise);
Z_oe = [y_oe, x];
m3_oe = oefit(Z_oe, nfit);
m3_arx = arxfit(Z_oe, nfit);
fprintf('n = %d\n', n)
disp('ARX expect: [0.2, -0.3, 0.4, -0.2]')
disp(m3_arx.theta')  % This should be the expected parameters defining yf
disp('OE expect: [0.2, -0.3, 0.4, -0.2]')
disp(m3_oe.theta')  % This should be the expected parameters defining yf

%% 3
clear all; close all;
disp('3')
% load('exercise1.mat'); z1 = [y(1:end/2,1), u(1:end/2, 1)]; z2 = [y(end/2+1:end,1), u(end/2+1:end, 1)];
load('exercise2.mat')
n_max = 12;

result_arx = analyze_model(@arxfit, z1, z2, n_max);
result_oe = analyze_model(@oefit, z1, z2, n_max);

H_arx = id2tf(result_arx.m);
H_oe = id2tf(result_oe.m); 

% Let plotting begin
idcompare(result_arx.m, z2, 1)
idcompare(result_oe.m , z2, 1)

figure; % ARX
subplot(2,2,1); step(H_arx); %step response
subplot(2,2,2); impulse(H_arx); % impulse response
subplot(2,2,3); pzmap(H_arx); % impulse response
subplot(2,2,4); nyquist(H_arx); % impulse response

figure; % OE
subplot(2,2,1); step(H_oe); %step response
subplot(2,2,2); impulse(H_oe); % impulse response
subplot(2,2,3); pzmap(H_oe); % impulse response
subplot(2,2,4); nyquist(H_oe); % impulse response

%% functions

function result = analyze_model(mfit, zt, zv, n_max)
    result = {};
    n1 = 1:n_max; n2 = 1:n_max; n3 = 0:n_max;
    n = {n1, n2, n3};
    [models, mq] = model_qualities(mfit, n, zt, zv);
    
    [~, I] = min(mq(:));
    [I1, I2, I3] = ind2sub(size(mq), I);
    
    n_opt = [n{1}(I1), n{2}(I2), n{3}(I3)];
    m_opt = mfit(zt, n_opt);
    
    result.index = [I1, I2, I3];
    result.models = models;
    result.model_qualities = mq;
    result.m = m_opt;
    result.n = n_opt;
end

function [models, mq] = model_qualities(mfit, n, z_train, z_valid)
    n1 = n{1}; d1 = length(n1);
    n2 = n{2}; d2 = length(n2);
    n3 = n{3}; d3 = length(n3);
    mq = zeros(d1, d2, d3);
    models = cell(d1, d2, d3);
    for i = 1:d1
        for j = 1:d2
            for k = 1:d3
                m = mfit(z_train, [n1(i), n2(j), n3(k)]);
                models{i,j,k} = m;
                y = z_valid(:, 1);
%                 y_sim = idsimulate(m, z_valid);
                y_sim = idpredict(m, z_valid, 1);
                mq(i,j,k) = immse(y, y_sim);
            end
        end
    end
end