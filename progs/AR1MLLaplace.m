% =========================================================================
% Estimate Laplace AR(1) model with constant and known variance of errors
% with Maximum Likelihood on simulated data
% =========================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =========================================================================
clearvars; clc; close all;
y = xlsread('../data/Laplace.xlsx',1); % load data
p = 1;                                 % set number of lags
const = 1;                             % model with constant
alph = 0.05;                           % significance level
MLLaPlace = ARpMLLaplace(y,p,const,alph); % estimate model using ARpMLLaplace

% Display results and compare to true values
TrueVals = [1; 0.8];
result = table(MLLaPlace.thetatilde,TrueVals);
result.Properties.VariableNames = {'ML_Estimate','True_Values'};
result.Properties.RowNames = {'c','\phi'};
disp(result)

% Compare to Gaussian ML estimates
MLGaussian = ARpMLLaplace(y,p,const,alph); % estimate model using ARpML
disp([MLLaPlace.thetatilde MLGaussian.thetatilde TrueVals]);