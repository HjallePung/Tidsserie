clear all;
close all;
projectData = load("projectData26.mat");
time = projectData.timeData;
waterflow = projectData.waterflowData;
waterflowData = iddata(waterflow);

lags = 360;

aaACFPACFN(waterflowData.OutputData, lags, "waterflow", 0.05)
