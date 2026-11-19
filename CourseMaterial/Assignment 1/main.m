%% 2.1 initial setup
clc;
clear;
close all;

A1 = [1, -1.79,  0.84]; %Polynomials
C1 = [1, -0.18, -0.11];

A2 = [1, -1.79];
C2 = [1, -0.18, -0.11];

arma1 = idpoly(A1,[],C1); %Create models
arma2 = idpoly(A2,[],C2);

rng(0)                     %Constants
sigma2 = 1.5;
m=100;
N = 300;
e = sqrt(sigma2) * randn(N,1);

y1 = filter(arma1.c,arma1.a,e); %Simulate data
y1 = y1(101:end);
y2 = filter(arma2.c,arma2.a,e);
y2 = y2(101:end);

r_theo = kovarians(arma1.c,arma1.a ,m);  %Plot covariance functions
stem(0:m,r_theo*sigma2)
hold on
r_est =covf(y1 ,m+1);
stem(0:m, r_est,'r')
legend('Theoretical','Estimated')
title('Covariance functions')

%% Plots
close all;      %Poles and zeroes
figure
subplot(121)
pzmap(arma1)
subplot(122)
pzmap(arma2)

figure
subplot(211)
plot(y1)
subplot(212)
plot(y2)

%% Estimate model
close all;
aaACFPACFN(y1,m,'y1',0.05)

%% Create model
clc;
na = 1;

data = iddata(y1);

ar_model = arx(data, na);

present(ar_model);

e_hat = aaFilter(ar_model.a ,1,data.OutputData,na);

plot(e_hat);

%% Re-estimate model
close all;
aaACFPACFN(e_hat,m,'residual',0.05)

%% 
nc = 1;

arma_model = armax(data,[na nc]);

present(arma_model);

e_hat = aaFilter(arma_model.a,arma_model.c, data.OutputData, na);

plot(e_hat);

%% Re-estimate
close all;
aaACFPACFN(e_hat,m,'residual',0.05);

%%
na = 2;

arma_model = armax(data,[na nc]);

present(arma_model);

e_hat = aaFilter(arma_model.a,arma_model.c, data.OutputData, na);

plot(e_hat);


%% Re-estimate
close all;
aaACFPACFN(e_hat,m,'residual',0.05);

%%
nc = 2;

arma_model = armax(data,[na nc]);

present(arma_model);

e_hat = aaFilter(arma_model.a,arma_model.c, data.OutputData, na);

plot(e_hat)

%% 2.2
clc;
clear;
close all;
noise = load("Code\data\noise.dat");
data_og = load("Code\data\data.dat");
data = iddata(data_og);
ar1 = arx(data,1);      %Models
ar2 = arx(data,2);
ar3 = arx(data,3);
ar4 = arx(data,4);
ar5 = arx(data,5);

rar1 = resid(data,ar1);     %Residual
rar1 = rar1(1+1:end);
rar2 = resid(data,ar2);
rar2 = rar2(2+1:end);
rar3 = resid(data,ar3);
rar3 = rar3(3+1:end);
rar4 = resid(data,ar4);
rar4 = rar4(4+1:end);
rar5 = resid(data,ar5);
rar5 = rar5(5+1:end);

figure;
subplot(231)
resid(data,ar1);
subplot(232)
resid(data,ar2);
subplot(233)
resid(data,ar3);
subplot(234)
resid(data,ar4);
subplot(235)
resid(data,ar5);

figure
subplot(231)                %Noise and Residual
plot(rar1.y)
hold on;
plot(noise(1+1:end))
legend("Residual", "Noise")
title("Noise and Residual (1)")

subplot(232)
plot(rar2.y)
hold on;
plot(noise(2+1:end))
legend("Residual", "Noise")
title("Noise and Residual (2)")

subplot(233)
plot(rar3.y)
hold on;
plot(noise(3+1:end))
legend("Residual", "Noise")
title("Noise and Residual (3)")

subplot(234)
plot(rar4.y)
hold on;
plot(noise(4+1:end))
legend("Residual", "Noise")
title("Noise and Residual (4)")

subplot(235)
plot(rar5.y)
hold on;
plot(noise(5+1:end))
legend("Residual", "Noise")
title("Noise and Residual (5)")

aaACFPACFN(rar1.y,50,"ar(1)",0.05); %Plot Acf and Pacf
aaACFPACFN(rar2.y,50,"ar(2)",0.05);
aaACFPACFN(rar3.y,50,"ar(3)",0.05);
aaACFPACFN(rar4.y,50,"ar(4)",0.05);
aaACFPACFN(rar5.y,50,"ar(5)",0.05);

% present(ar1);
% present(ar2);
% present(ar3);
% present(ar4);
% present(ar5);

%% ARMA
clc;
close all;
am11 = armax(data,[1 1]);  %ARMA models
am12 = armax(data,[1 2]);
am21 = armax(data,[2 1]);
am22 = armax(data,[2 2]);

res11 = resid(data,am11);   %Residuals
res11 = res11(1+1:end);
res12 = resid(data,am11);
res12 = res12(1+1:end);
res21 = resid(data,am11);
res21 = res21(2+1:end);
res22 = resid(data,am11);
res22 = res22(2+1:end);

subplot(221)    %Plot noise and residuals
plot(res11.y)
hold on;
plot(noise(1+1:end))
legend("Residual", "Noise")
title("Noise and Residual (11)")

subplot(222)
plot(res12.y)
hold on;
plot(noise(1+1:end))
legend("Residual", "Noise")
title("Noise and Residual (12)")

subplot(223)
plot(res21.y)
hold on;
plot(noise(2+1:end))
legend("Residual", "Noise")
title("Noise and Residual (21)")

subplot(224)
plot(res22.y)
hold on;
plot(noise(2+1:end))
legend("Residual", "Noise")
title("Noise and Residual (22)")

aaACFPACFN(res11.y,50,"arma(1,1)",0.05);
aaACFPACFN(res12.y,50,"arma(1,2)",0.05);
aaACFPACFN(res21.y,50,"arma(2,1)",0.05);
aaACFPACFN(res22.y,50,"arma(2,2)",0.05);

% present(am11);
% present(am21);
% present(am12);
% present(am22);

%% 2.3
clc;
clear;
close all;

rng(0)              %Simulate Data
A = [1 -1.5 0.7];
C = [1 zeros(1,11) -0.5];
A12 = [1 zeros(1,11) -1];
Astar = conv(A, A12);
e = randn(600, 1);
y = filter(C, Astar, e);
y = y(101:end);

plot(y)             
xlabel('Time')
ylabel('y')
title('Simulated Output')
grid on

aaACFPACFN(y,50,"data",0.05);

%% Remove seasonality
ys = filter(A12, 1, y);
ys = ys(length(A12):end);
data = iddata(ys);

%Ar(2)
modelinit = idpoly([1 0 0], [], []);
modelarmax = pem(data, modelinit);

e_hat = aaFilter(modelarmax.a,modelarmax.c,data.OutputData,20);

aaACFPACFN(e_hat, 50,"Residual",0.05);

%% Add MA part
clc;
modelinit = idpoly([1 0 0], [], [1 zeros(1,12)]);
modelinit.Structure.c.Free = [zeros(1,12) 1];
modelarmax = pem(data, modelinit);

present(modelarmax)

e_hat = aaFilter(modelarmax.a,modelarmax.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);
%%
clc;
close all;

rng(0);             %With 10000 samples
A = [1 -1.5 0.7];
C = [1 zeros(1,11) -0.5];
A12 = [1 zeros(1,11) -1];
Astar = conv(A, A12);
e = randn(10000, 1);
y = filter(C, Astar, e);
y = y(101:end);

% plot(y)
% xlabel('Time')
% ylabel('y')
% title('Simulated Output')
% grid on

% aaACFPACFN(y,20,"data",0.05);

ys = filter(A12, 1, y);
ys = ys(length(A12):end);

data = iddata(ys);

modelinit = idpoly([1 0 0], [], []);

modelarmax = pem(data, modelinit);

e_hat = aaFilter(modelarmax.a,modelarmax.c,data.OutputData,20);

% aaACFPACFN(e_hat, 20,"Residual",0.05)

modelinit = idpoly([1 0 0], [], [1 zeros(1,12)]);
modelinit.Structure.c.Free = [zeros(1,12) 1];

modelarmax = pem(data, modelinit);

present(modelarmax)

e_hat = aaFilter(modelarmax.a,modelarmax.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

%% Without removing season
clc;
close all;
rng(0)

A = [1 -1.5 0.7];
C = [1 zeros(1,11) -0.5];
A12 = [1 zeros(1,11) -1];
Astar = conv(A, A12);
e = randn(600, 1);
y = filter(C, Astar, e);
y = y(101:end);

data2 = iddata(y);

% modelinit = idpoly([1 zeros(1,11) 1], [], []);
% modelinit.Structure.a.Free =[1 1 zeros(1,10) 1];

modelinit = idpoly([1 0], [], []);

mod2 = pem(data2, modelinit);

e_hat = aaFilter(mod2.a,mod2.c,data2.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2);

%%

modelinit2 = idpoly([1 0 0], [], []);
% modelinit2.Structure.a.Free = [1 1 1 zeros(1,9) 1];

mod2 = pem(data2, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data2.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2)

%%

modelinit2 = idpoly([1 0 0], [], [1 zeros(1,12)]);
modelinit2.Structure.c.Free = [1 zeros(1,11) 1];

mod2 = pem(data2, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data2.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2)


%%
clc;
close all;
modelinit2 = idpoly([1 zeros(1,14)], [], [1 zeros(1,12)]);

modelinit2.Structure.a.Free = [1 1 1 zeros(1,9) 1 1 1];
modelinit2.Structure.c.Free = [1 zeros(1,11) 1];

mod2 = pem(data2, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data2.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

figure
plot(resid(mod2,data2).y)
hold on;
plot(resid(modelarmax,data).y)
legend("Seasonality not removed" ,"Seasonality removed")
title("Residuals")

present(mod2)

%% Svedala
clc;
clear;
close all;
y = load("Code\data\svedala.mat").svedala;
data = iddata(y);

figure
plot(data)

aaACFPACFN(data.OutputData,50,"svedala",0.05);

%%
modelinit3 = idpoly([1 0 0], [], []);

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,50);

present(mod3);

aaACFPACFN(e_hat,50,"residual",0.05);

%%
modelinit3 = idpoly([1 zeros(1,24)], [], [1 zeros(1,23)]);
modelinit3.Structure.a.Free = [1 1 1 zeros(1,7) 1 1 zeros(1,11) 1 1];
modelinit3.Structure.c.Free = [1 zeros(1,22) 1];
mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,50);

figure;
plot(e_hat)

present(mod3);

aaACFPACFN(e_hat,250,"residual",0.05);

%% Svedala removing seasonality
A24 =[1 zeros(1,25) -1];
ys = aaFilter(A24,1,y,50);

data = iddata(ys);

aaACFPACFN(data.OutputData, 50,"no season", 0.05)

modelinit3 = idpoly([1 0 0], [], []);

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,50);

aaACFPACFN(e_hat, 50,"no season", 0.05)

present(mod3);

%%
modelinit3 = idpoly([1 0 0], [], [1 zeros(1,23) 1]);
modelinit3.Structure.c.Free = [1 zeros(1,23) 1];

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,50);

aaACFPACFN(e_hat, 50,"no season", 0.05)

present(mod3);

%%
modelinit3 = idpoly([1 zeros(1,24)], [], [1 zeros(1,23) 1]);
modelinit3.Structure.a.Free = [1 1 1 zeros(1,21) 1 ];
modelinit3.Structure.c.Free = [1 zeros(1,23) 1];

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,50);

aaACFPACFN(e_hat, 50,"no season", 0.05)

present(mod3);
