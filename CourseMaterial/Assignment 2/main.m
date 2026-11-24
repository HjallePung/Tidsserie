clear all;
close all;
rng(0)
n = 500; % Number of samples

A3 = [1 .5];
C3 = [1 -.3 .2];
w = sqrt(2) * randn(n+100, 1);
x = filter(C3, A3, w); % Create the input

A1 = [1 -.65];
A2 = [1 .90 .78];
C = 1;
B = [0 0 0 0 .4];

e = sqrt(1.5) * randn(n+100, 1);
y = filter(C, A1, e) + filter(B, A2, x); % Create the output

x = x(101:end);, y = y(101:end); % Omit initial samples

clear A1;, A2;, C;, B;, e;, w;, A3;, C3;

modelinit = idpoly([1 0], [], []);

mod = pem(x, modelinit);

e_hat = aaFilter(mod.a,mod.c,x,50);

aaACFPACFN(e_hat, 50,"x", 0.05);

present(mod);


%%
epst = aaFilter(mod.a, mod.c, y, 50);
wt = aaFilter(mod.a, mod.c, x, 50);

%%
aaCrossCorr(wt, epst, n, 40)
%%
A2 = [1 0 0];
B = [0 0 0 0 1];

Mi = idpoly([1], [B], [], [], [A2]);
z = iddata(y, x);

Mba2 = pem(z, Mi);
present(Mba2)

etilde = resid(Mba2, z);

aaACFPACFN(etilde.OutputData, 40, "etilde", 0.05)

aaCrossCorr(x, etilde.OutputData, n, 40)

%%
xfilt = filter(Mba2.b, Mba2.f, x);
figure;
plot([y(length(A2):end) xfilt(length(A2):end)])
legend('y', 'x_{filtered}')

%%
% Subtract and plot to check that we did correct
figure;
plot([y(length(A2):end) y(length(A2):end)-xfilt(length(A2):end)])
legend('y', 'subtracted values')

%%
% Subtract the model based on x from y
etTilde = y(length(A2):end)-xfilt(length(A2):end);

aaACFPACFN(etTilde, 40, "etTilde", 0.05)

%%
modelinit1 = idpoly([1 0], [], []);

mod2 = pem(etTilde, modelinit1);

present(mod2)

e = resid(mod2, etTilde);

aaACFPACFN(e, 40, "AR(1)", 0.05)
%%
% Check crosscovariance between etTilde and x
aaCrossCorr(x, e, 498, 40)
%%
A1 = [1 0];
A2 = [1 0 0];
B = [0 0 0 0 1];
C = [1];

Mi = idpoly(1, B, C, A1, A2);
z = iddata(y, x);

MboxJ = pem(z, Mi);
present(MboxJ)

ehat = resid(MboxJ, z);

aaACFPACFN(ehat.OutputData, 40, "full model", 0.05)
%%
aaCrossCorr(ehat.OutputData, x, 498, 40)
%%
xfilt = filter(Mi.b, Mi.f, x);
figure;
plot([y(length(A2):end) xfilt(length(A2):end)])
legend('y', 'x_{filtered}')
%%
xfilt = filter(Mi.b, Mi.f, x);
figure;
plot([y(length(A2):end) y(length(A2):end)-xfilt(length(A2):end)])
legend('y', 'subtracted')

%%
clear all;
close all;
y = load("CourseMaterial\Code\data\svedala.mat").svedala;
fprintf('Mean of y: %.4f\n', mean(y));
y = y - mean(y);

plot(y)

aaACFPACFN(y, 40, "y", 0.05)

A = [1 -1.79 0.84];
C = [1 -0.18 -0.11];

k = 3;
[Fk3 Gk3] = polydivision(C, A, k);

yhat_3 = aaFilt(Gk3,C,y);

error_3 = y(k+1:end) - yhat_3;

k = 26;
[Fk26 Gk26] = polydivision(C, A, k);

yhat_26 = aaFilt(Gk26,C,y);

error_26 = y(4:end) - yhat_26;

aaACFPACFN(error_3, 40, "error", 0.05)
%%
fprintf('Mean error of k = 3: %.4f\n', mean(error_3));
fprintf('Mean error of k = 26: %.4f\n', mean(error_26));


%%
sigma3_2_e = var(aaFilt(A, C, y));
sigma26_2_e = var(aaFilt(A, C, y));

var_theoretical3 = sigma3_2_e * sum(Fk3.^2);
var_theoretical26 = sigma26_2_e * sum(Fk26.^2);

var_estimated3 = var(error_3);
disp('for 3:')
fprintf('Theoretical variance: %.4f\n', var_theoretical3);
fprintf('Estimated variance:   %.4f\n', var_estimated3);

fprintf('Relative difference: %.2f%%\n', ...
    100 * (var_estimated3 - var_theoretical3) / var_theoretical3);

disp('for 26:')
var_estimated26 = var(error_26);
fprintf('Theoretical variance: %.4f\n', var_theoretical26);
fprintf('Estimated variance:   %.4f\n', var_estimated26);

fprintf('Relative difference: %.2f%%\n', ...
    100 * (var_estimated26 - var_theoretical26) / var_theoretical26);

%%
c3 = 1.96*sqrt(sigma3_2_e*sum(Fk3.^2));
c26 = 1.96*sqrt(sigma26_2_e*sum(Fk26.^2));
fprintf('3 95%% confidence intervall: +-%.4f\n', c3)
fprintf('26 95%% confidence intervall: +-%.4f\n', c26)

sum(error_3>c3)
sum(error_26>c26)

%%
figure
plot(yhat_3(24:end))
hold on;
plot(y(27:end))
hold on;
plot(yhat_26)
legend('yhat_3', 'y', 'yhat_{26}');

figure
plot(error_3)
hold on;
plot(error_26)
legend('error_3', 'error_{26}');

%%
aaACFPACFN(error_3, 40, 'error_3', 0.05);
aaACFPACFN(error_26, 40, 'error_{26}', 0.05);

%%
close all;
clear all;
x = load("CourseMaterial\Code\data\sturup.mat").sturup;
y = load("CourseMaterial\Code\data\svedala.mat").svedala;

A = [1 -1.49 0.57];
B = [0 0 0 0.28 -0.26];
C = [1];


k = 3; %First case
[Fk3 Gk3] = polydivision(C, A, k);

BF3 = conv(Fk3, B);

[Fk3_hat Gk3_hat] = polydivision(BF3, C, k);

xhat_3 = filter(Gk3,C,x);

yhat_3 = aaFilter(Fk3_hat,1,x,26) + aaFilter(Gk3_hat,C,x,26) + aaFilter(Gk3,C,y, 26);

error_3 = y(27:end) - yhat_3;

k = 26;

[Fk26 Gk26] = polydivision(C, A, k);

[Fk26 Gk26] = polydivision(C, A, k);

BF26 = conv(Fk26, B);

[Fk26_hat Gk26_hat] = polydivision(BF26, C, k);

yhat_26 = aaFilter(Fk26_hat,1,x,26) + aaFilter(Gk26_hat,C,x,26) + aaFilter(Gk26,C,y, 26);

error_26 = y(27:end) - yhat_26;

aaACFPACFN(error_3, 40, "error_3", 0.05)
aaACFPACFN(error_26, 40, "error_{26}", 0.05)

%% plot the predictions
figure;
plot(y(27:end))
hold on;
plot(yhat_3)
%hold on;
%plot(yhat_26)
legend('y' ,'y_3')

%% Variance of y3
figure;
plot(error_3)

fprintf("Variance for k = 3: %f", var(error_3))

%% Without the Fhat term
y3_wrong = yhat_3 - aaFilter(Fk3_hat,1,x,26);

figure;
plot(y(27:end))
hold on;
plot(y3_wrong)

error_wrong = y(27:end) - y3_wrong;

figure;
plot(error_wrong)

fprintf("Variance for k = 3: %f", var(error_wrong))