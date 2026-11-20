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
A2 = [1 0 0 0 0];
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
legend('y', 'x_filtered')

%%
% Subtract and plot to check that we did correct
figure;
plot([y(length(A2):end) y(length(A2):end)-xfilt(length(A2):end)])
legend('y', 'subtracted values')

