%% 2.1 initial setup
clc;
clear;
close all;

A1 = [1, -1.79,  0.84];
C1 = [1, -0.18, -0.11];

A2 = [1, -1.79];
C2 = [1, -0.18, -0.11];

arma1 = idpoly(A1,[],C1);
arma2 = idpoly(A2,[],C2);

rng(0)
sigma2 = 1.5;
m=20;
N = 300;
e = sqrt(sigma2) * randn(N,1);

y1 = filter(arma1.c,arma1.a,e);
y1 = y1(101:end);
y2 = filter(arma2.c,arma2.a,e);
y2 = y2(101:end);

r_theo = kovarians(arma1.c,arma1.a ,m);
stem(0:m,r_theo*sigma2)
hold on
r_est =covf(y1 ,m+1);
stem(0:m, r_est,'r')
legend('Theoretical','Estimated')
title('Covariance functions')

%% Plots
close all;
figure
pzmap(arma1)
figure
pzmap(arma2)

figure
subplot(211)
plot(y1)
subplot(212)
plot(y2)

%% Re-estimate model
close all;
aaACFPACFN(y1,m,'y1',0.05)

%% 
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
close all;
noise = load("Code\data\noise.dat");
data_og = load("Code\data\data.dat");
data = iddata(data_og);
ar1 = arx(data,1);
ar2 = arx(data,2);
ar3 = arx(data,3);
ar4 = arx(data,4);
ar5 = arx(data,5);

rar1 = resid(data(1+1:end),ar1);
rar2 = resid(data(2+1:end),ar2);
rar3 = resid(data(3+1:end),ar3);
rar4 = resid(data(4+1:end),ar4);
rar5 = resid(data(5+1:end),ar5);

figure;
resid(data(5+1:end),ar5);
figure;
resid(data(4+1:end),ar4);

figure
plot(rar1.y)
hold on;
plot(noise(2:end))
legend("Residual", "Noise")

figure
plot(rar4.y)
hold on;
plot(noise(5:end))

present(ar1);
present(ar2);
present(ar3);
present(ar4);
present(ar5);

%% ARMA
clc;
close all;
am11 = armax(data,[1 1]);
am12 = armax(data,[1 2]);
am21 = armax(data,[2 1]);
am22 = armax(data,[2 2]);

resid(data(1+1:end),am11);

present(am11);
present(am21);
present(am12);
present(am22);

%% 2.3
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

plot(y)
xlabel('Time')
ylabel('y')
title('Simulated Output')
grid on

aaACFPACFN(y,20,"data",0.05);

%% 
ys = filter(A12, 1, y);

ys = ys(length(A12):end);

data = iddata(ys);

modelinit = idpoly([1 0 0], [], []);

modelarmax = pem(data, modelinit);

e_hat = aaFilter(modelarmax.a,modelarmax.c,data.OutputData,20);

aaACFPACFN(e_hat, 20,"Residual",0.05)

%% 
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

rng(0)

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

data = iddata(y);

modelinit = idpoly([1 zeros(1,11) 1], [], []);

modelinit.Structure.a.Free =[1 1 zeros(1,10) 1];

mod2 = pem(data, modelinit);

e_hat = aaFilter(mod2.a,mod2.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2);

%%

modelinit2 = idpoly([1 0.1 zeros(1,10) 0.1], [], []);

modelinit2.Structure.a.Free = [1 1 1 zeros(1,9) 1];

mod2 = pem(data, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2)

%%

modelinit2 = idpoly([1 0.1 zeros(1,10) 0.1], [], [1 1]);

modelinit2.Structure.a.Free = [1 1 1 zeros(1,9) 1];

mod2 = pem(data, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

present(mod2)


%%
clc;
close all;
modelinit2 = idpoly([1 0.1 0.1 zeros(1,9) 0.1], [], [1 1 1]);

modelinit2.Structure.a.Free = [1 1 1 1 zeros(1,8) 1];

mod2 = pem(data, modelinit2);

e_hat = aaFilter(mod2.a,mod2.c,data.OutputData,20);

aaACFPACFN(e_hat,20,"residual",0.05);

figure
plot(resid(mod2,data).y)
hold on;
plot(resid(modelarmax,ys).y)

present(mod2)

%% Svedala
clc;
clear;
close all;
y = load("Code\data\svedala.mat").svedala;
data = iddata(y);

figure
plot(data)

aaACFPACFN(data.OutputData,250,"svedala",0.05);

modelinit3 = idpoly([1 1 1], [], []);

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,250);

present(mod3);

aaACFPACFN(e_hat,250,"residual",0.05);

%% Svedala removing seasonality
A25 =[1 zeros(1,24) -1];
ys = aaFilter(A25,1,y,250);

data = iddata(ys);

aaACFPACFN(data.OutputData, 250,"no season", 0.05)

modelinit3 = idpoly([1 1 1], [], []);

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,250);

aaACFPACFN(e_hat, 250,"no season", 0.05)

present(mod3);

%%
aaACFPACFN(data.OutputData, 250,"no season", 0.05)

modelinit3 = idpoly([1 1 1], [], [1 zeros(1,24) 1]);

modelinit3.Structure.c.Free = [1 zeros(1,24) 1];

mod3 = pem(data, modelinit3);

e_hat = aaFilter(mod3.a,mod3.c,data.OutputData,250);

aaACFPACFN(e_hat, 250,"no season", 0.05)

present(mod3);
