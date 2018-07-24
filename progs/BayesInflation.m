% =======================================================================
% Bayesian estimation of an AR(2) model of quarterly inflation
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
% =======================================================================

clearvars; clc;close all;

%% DATA HANDLING
data = xlsread('../data/data.xlsx','QuarterlyInflation'); % Load data
T=size(data,1);                        % Determine sample length, i.e. how many quarters
X = [ones(T,1) lagmatrix(data,1:2)];   % Matrix of regressors, i.e. c, y(t-1) and y(t-2) for AR(2) model with constant
X=X(3:end,:);                          % Remove the first 2 observations
y=data(3:end,1);                       % Remove the first 2 observations of dependent variable
T=size(y,1);                           % Sample length after adjustment

%% Set priors
% Priors for B (normal distribution):
B0=zeros(3,1);   % Prior mean for coefficients
Sigma0=eye(3);   % Prior variance for coefficients
%Priors for inv_sigma2=1/sigma^2 (inverse Gamma distribution):
s0=1;            % Prior degrees of freedom equal to number of variables
v0=0.1;          % Prior scale parameter

%% Options for Gibbs sampler
nreps=5000;      % Total number of Gibbs iterations
nburn=4000;      % Number of burn-in iterations

%% Initialization of Gibbs sampler
inv_sigma2=1;               % Initialize 1/sigma^2
out1=zeros(nreps-nburn,3);  % Coefficient estimates
out2=zeros(nreps-nburn,1);  % Precision estimates

%% Gibbs sampling
count=1; 
for jj=1:nreps
    % Sample B conditional on sigma^2 from N(M,V)
    invSigma0 = Sigma0\eye(size(Sigma0,1));
    % Posterior mean of B
    M=(invSigma0+(1/inv_sigma2)*(X'*X))\(invSigma0*B0+(1/inv_sigma2)*X'*y);
    % Posterior variance of B
    V=inv(invSigma0+(1/inv_sigma2)*(X'*X));

    chck=-1; %Check for stability of draw (to avoid explosive draws)
    while chck<0
        B=mvnrnd(M,V)';      % Take a draw from multivariate normale
        A=[B(2) B(3); 1 0];  % Companion matrix
        ee=max(abs(eig(A))); 
        if ee<=1
            chck=1;          % AR model is stable if the eigenvalues are less than or equal to 1 in absolute value
        end
    end

    % Sample sigma2 conditional on B from IG(s1,v1)
    resids=y-X*B;           %Compute residuals conditional on B
    s1=s0+T;                %Compute posterior degrees of freedom 
    v1=v0+resids'*resids;   %Compute posterior scale matrix
    inv_sigma2 = 1/gamrnd(s1,1/v1); %draw from inverse gamma

    % Save draws for inference if burn-in phase is passed
    if jj>nburn
        out1(count,:)= B;
        out2(count,:)= inv_sigma2;
        count=count+1;
    end
end

%% Plot marginal posterior distributions
x = -1:0.1:1; % values to plot normal pior distribution
xx = 2:0.1:5; % values to plot gamma pior distribution
figure('name','Marginal Posterior Distributions');

subplot(2,2,1)
histogram(out1(:,1),50,'Normalization','pdf');
axis tight
hold on;
plot(x,normpdf(x,B0(1),sqrt(Sigma0(1,1)))); % Normal prior distribution
title('Constant')
hold off;

subplot(2,2,2)
histogram(out1(:,2),50,'Normalization','pdf');
axis tight
hold on;
plot(x,normpdf(x,B0(2),sqrt(Sigma0(2,2)))); % Normal prior distribution
title('AR(1) coefficient')
hold off;

subplot(2,2,3)
histogram(out1(:,3),50,'Normalization','pdf');
axis tight
hold on;
plot(x,normpdf(x,B0(3),sqrt(Sigma0(3,3)))); % Normal prior distribution
title('AR(2) coefficient')
hold off;

subplot(2,2,4)
histogram(1./out2(:,1),50,'Normalization','pdf');
axis tight
hold on;
plot(xx,gampdf(xx,s0,1/v0)); % Gamma prior distribution
title('\sigma^2')
hold off;

%% Compute and display estimates from posterior
MB1=mean(out1);              % Compute mean of the marginal posterior distribution of B
SB1=std(out1);               % Compute standard error of B
EB1=prctile(out1,[5 95]);    % Compute 95% error bands of B

MB2=mean(out2);              % Compute mean of the marginal posterior distribution of 1/sigma^2
SB2=std(out2);               % Compute standard error of 1/sigma^2
EB2=prctile(out2,[5 95]);    % Compute 95% error bands of 1/sigma^2

display(['Mean of B:           ' num2str(MB1)])
display(['Standard error of B: ' num2str(SB1)])
display(['OLS estimates:       ' num2str(((X'*X)\(X'*y))')])