% =========================================================================
% Example on how to load and visualize data with plot, histfit, normplot
% =========================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =========================================================================
clearvars; clc;close all;
%% Load data
% load data from different folders, e.g. '../' goes one subdirectory down
data_raw=xlsread('../data/Quarterly_GDP_Norway.xlsx',1);

dates=data_raw(:,1);
data=data_raw(:,2);

%% Figures
f1=figure('name','Data');              % open new window for the plots
h1=plot(dates,data,'linewidth',2);     % simple plot

f2=figure('name','Histogram of data'); % open new window for the plots
h2=histfit(data,10,'normal');          % histogram, fitted normal distribution

f3=figure('name','Normplot of data');  % open new window for the plots
h3=normplot(data);                     % Q-Q plot

%% Estimates
mean(data) % empirical estimation of average and standard deviation
std(data)  % empirical estimation of standard deviation