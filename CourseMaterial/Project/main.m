%% Initialize Data
clc;
clear all;
close all;
projectData = load("projectData26.mat");
time = projectData.timeData;
waterflow = projectData.waterflowData;
waterflowData = iddata(waterflow);
%%
lags = 20;

plot(waterflow)

aaACFPACFN(waterflowData.OutputData, lags, "waterflow", 0.05);

%%  Plot data

waterflowDataZeroMean = iddata(waterflowData.OutputData - mean(waterflowData.OutputData));

aaACFPACFN(waterflowDataZeroMean.OutputData, lags, "waterflow", 0.05);

%% Extract 6 weeks
waterflowTrainingData =waterflowDataZeroMean(500: 500 + 41);

aaACFPACFN(waterflowTrainingData.OutputData, lags, "waterflow", 0.05);

%% Create Model
modelinit = idpoly([1 0], [], []);
% modelinit3.Structure.a.Free = [1 1 1 zeros(1,7) 1 1 zeros(1,11) 1 1];
mod1 = pem(waterflowTrainingData, modelinit);

present(mod1)

e_hat = aaFilter(mod1.a,mod1.c,waterflowTrainingData.OutputData,lags);

aaACFPACFN(e_hat, lags, "residual", 0.05);

figure
plot(e_hat)

%% Model Step
modelinit = idpoly([1 0], [], [1 0]);
% modelinit.Structure.c.Free = [1 zeros(1,4) 1];
mod1 = pem(waterflowTrainingData, modelinit);

present(mod1)

e_hat = aaFilter(mod1.a,mod1.c,waterflowTrainingData.OutputData,lags);

aaACFPACFN(e_hat, lags, "residual", 0.05);

figure
plot(e_hat)

