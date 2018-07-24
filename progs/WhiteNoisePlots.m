% =========================================================================
% Generate and plot white noise and moving-average processes
% =========================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =========================================================================
clearvars; clc;close all;
%% Generate white noise
sigma=1;            % set value for standard deviation
T=200;              % set value for number of observations
epsi = randn(T,1);  % draw Tx1 vector of Gaussian random variables
size(epsi)          % check size of epsi, should be Tx1

y1 = nan(T,1);      % initialize a Tx1 vector with nan (not a number)
for t=1:T           % the variable t runs from 1,2....,T
    y1(t,1) = epsi(t)*sigma;
    % epsi are drawn from N(0,1), we can scale the standard error 
    % by multiplying with sigma
end

%% Generate moving average
windowSize=5;       % set value for window of moving average
y2 = nan(T,1);      % initialize output vector with nan
for t=3:T-2         % since epsi is Tx1, t cannot start at 1 as we need epsi(t-2);
                    % t cannot end at T, as we end the loop with epsi(t+2);
    y2(t) = 1/windowSize*(epsi(t-2)+epsi(t-1)+epsi(t)+epsi(t+1)+epsi(t+2));
end

%% Plots
f=figure('name','White Noise Plots');   % open new window for figure
subplot(1,2,1);                         % subplot(rows,columns, plotindex)
h1=plot(y1);                            % plot white noise
title('White Noise');                   % set title
subplot(1,2,2);                         % subplot(rows,columns, plotindex)
h2=plot(y2);                            % plot moving-average
title('5-point Moving Average');        % set title