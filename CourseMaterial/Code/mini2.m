clear; clc;
close all;
N = 100;

et = randn([N,1]);

figure
plot(et)

roots([1, - 0.79 , 0.84])

proc = filter([1, - 0.79, 0.84],[1, -2.81, -0.11],et);

figure
plot(proc)

figure
acf(proc,20, 0.05,1);

figure
pacf(proc,20, 0.05,1);

