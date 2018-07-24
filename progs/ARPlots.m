% =======================================================================
% Simulate and plot autoregressive processes and random walks
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================
clearvars; clc;close all;
%% Generate and plot autoregressive processes
phi=[-0.8 0.4 0.9 1.01];    % different values for the phi coefficient
sigma=1;                    % value for standard deviation of white noise
T=200;                      % value for number of observations
y=nan(T,4);                 % initiailize output vector with nan
y(1,:)=zeros(1,4);          % set first period to zero
for t=2:T
    for j=1:4
        y(t,j)=phi(j)*y(t-1,j)+randn()*sigma; % Simulate time series
    end
end

nameAr={'ARmin08','AR04','AR09','AR101'}; % create cell with titles of plots
                                          % note the use of curly brakets 
                                          % to deal with strings of
                                          % different length
figure('name','AR Plots');                % open new window for figure
for i=1:4
    subplot(2,2,i);
    plot(y(:,i));
    title(nameAr{i}); % use curly brackets to access cells with strings
end


%% Generate and plot random walks
nRW = 16;               % number of Random Walks to generate
sigma=1;                % value of standard error of white noise
T=200;                  % number of observations
y=nan(T,nRW);           % initialize output vector with nan
y(1,:)=zeros(1,nRW);    % set first period to zero
for t=2:T
    for j = 1:nRW
        y(t,j)=y(t-1,j)+randn()*sigma;
    end
end

figure('name','Random Walks');
for j=1:nRW
    subplot(4,4,j);
    plot(y(:,j));    
end