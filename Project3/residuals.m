fname = fullfile('data','residuals.mat');
load(fname);
Ts = 0.001;
samplePeriod = 100;
zt = iddata(YT(1:samplePeriod:end,1), UT(1:samplePeriod:end), samplePeriod*Ts);  % training data
zv = iddata(YV(1:samplePeriod:end,1), UV(1:samplePeriod:end), samplePeriod*Ts);  % validation data

V = ivstruc(zt,zv,struc(1:10,1:10,1:10));
nn = selstruc(V,0);


disp('Nonlinear Identification')
NL = wavenet('NumberOfUnits',5);
wavelet = nlarx(zt,nn,NL);
NL = sigmoidnet('NumberOfUnits',5);
sigmoidnet = nlarx(zt,nn,NL);
NL = treepartition('NumberOfUnits',5);
treepartition = nlarx(zt,nn,NL);
linearARX = arx(zt,nn);
linearModel = nlarx(zt,linearARX);
figure
compare(zv, linearARX, wavelet,...
    sigmoidnet, treepartition,linearModel, 5);
figure
compare(zv, linearARX, wavelet,...
    sigmoidnet, treepartition,linearModel, Inf);