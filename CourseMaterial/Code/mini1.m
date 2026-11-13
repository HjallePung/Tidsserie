clc;
clear;
close all;

[y,Fs] = audioread("data\fa.wav");

% sound(y)

figure
plot(y)

vowel = y(2400:2600);

figure
plot(vowel)

rho = acf(vowel,200);

figure
plot(rho)

%I think the main periodicity is something like 15Hz as this is is the
%distance between peaks in the acf graph
%Nah ok it seems to be 45 Hz, which is like the distance between the two
%largest peaks in the acf

%Power density function FFT of ACF
psd = abs(fft(rho)); 
figure
plot(psd)

%Mby different

psd = 1/length(vowel) * abs(fft(vowel)).^2;
figure
plot(psd)

x = linspace(1,9,1000);
figure
plot(cos(19*x) + 0.8*cos(28*x))
