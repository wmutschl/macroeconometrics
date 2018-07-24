% =======================================================================
% Bayesian Estimation of a VAR(2) for US data on Federal Funds Rate, 
% Government Bond Yield, Unemployment and Inflation from 2007m1 2010m12.
% Bayesian estimation with Gibbs Sampling using a Minnesota Prior for the
% VAR coefficients adjusted for the Zero-Lower-Bound and an inverse Wishart
% Prior for the covariance matrix.
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
% =======================================================================
clearvars; clc;close all;

%% PRELIMINARIES
% Define specification of the VAR model
const = 1;              % 0: no constant, 1: constant, 2: constant and linear trend
p = 2;                  % Number of lags on dependent variables

% Define Minnesota prior for BVAR model
prior_adjust   = 1;     % Set to 1 to adjust prior for Zero-Lower-Bound, 0 usual Minnesota prior
hyperparams(1) = 0.5;   % tightness parameter for Minnesota prior on lags of own variable
hyperparams(2) = 0.5;   % tightness parameter for Minnesota prior on lags of other variable
hyperparams(3) = 1;     % tightness of prior on exogenous variables, i.e. constant term, trends, etc

% Gibbs-related preliminaries
nsave = 10000;          % Final number of draws to keep
nburn = 30000;          % Draws to discard (burn-in)
ntot = nsave+nburn;     % Total number of draws
it_print = 1000;        % Print on the screen every "it_print"-th iteration

%% DATA HANDLING
% Load monthly US data on FFR, govt bond yield, unemployment and inflation
Yraw = xlsread('../data/data.xlsx','USZLB'); %'Yraw' is a matrix with T rows by K columns 
[Traw,K] = size(Yraw);                       % Get initial dimensions of dependent variable
Ylag = lagmatrix(Yraw,1:p);                  % Generate lagged Y matrix. This will be part of the Z matrix
% Now define matrix Z which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies), also get rid of NA observations
if const == 0
    Z = transpose(Ylag(p+1:Traw,:));  
elseif const == 1
    Z = transpose([ones(Traw-p,1) Ylag(p+1:Traw,:)]);
elseif const == 2
    Z = transpose([ones(Traw-p,1) transpose((p+1):Traw) Ylag(p+1:Traw,:)]);
end
Y = transpose(Yraw(p+1:Traw,:)); % Dependent variable in i-th equation, get rid of NA observations
[totcoeff,T] = size(Z);          % Get size of final matrix Z
ZZt = Z*Z';                      % auxiliary matrix product

%% PRIOR SPECIFICATION
% Get standard specification of Minnesota Normal-Inverse-Wishard Prior
[alpha_prior, V_prior, inv_V_prior, v_prior, S_prior, inv_S_prior] = BVARMinnesotaPrior(Yraw,const,p,hyperparams);
if prior_adjust % Manually adjust for Zero-Lower-Bound
    I_K = eye(K);
    alpha_prior(const*K+1: const*K+K*K) = 0.95*I_K(:); % AR prior instead of random walk
    % Adjust the prior variance of var coefficients to reflect ZLB
    tmp = zeros(K,K);
    tmp(1,2:K) = 1; %Focus on A1_12, A1_13, A1_14, A2_12, A2_13, A2_14
    Atmp = [zeros(K,1) tmp tmp];
    idx = find(Atmp==1);
    for j = idx
        V_prior(j,j) = 1e-9; % small for coefficients we want close to zero due to the ZLB
    end
    inv_V_prior = diag(1./diag(V_prior));
end

%% INITIALIZATION
A_OLS = (Y*Z')/ZZt;                                 % Get OLS estimators
alpha_OLS = A_OLS(:);                               % Vectorize
resid_OLS = Y - A_OLS*Z;                            % Compute OLS residuals
SIGMA_OLS = (resid_OLS*resid_OLS')./(T-K*p-const);  % OLS estimator for covariance matrix

% Initialize Bayesian posterior parameters using OLS values
alpha = alpha_OLS; % This is the first draw from the posterior of alpha
A = A_OLS;         % This is the first draw from the posterior of A
SIGMA = SIGMA_OLS; % This is the single draw from the posterior of SIGMA

%% START GIBBS SAMPLING
A_draws = zeros(nsave,K,totcoeff); % Storage space for posterior draws
SIGMA_draws = zeros(nsave,K,K);    % Storage space for posterior draws
tic;                               % Start timer
waitb = waitbar(0,'Number of iterations');
for irep = 1:ntot
    if mod(irep,it_print) == 0
        waitbar(irep/ntot); % Update waitbar every "it_print"-th step
    end
    %Posterior of alpha|SIGMA,Y ~ N(alpha_post,V_post)
    invSIGMA = inv(SIGMA);    
    V_post = (inv_V_prior+kron(ZZt,invSIGMA))\eye(size(inv_V_prior,1));
    alpha_post = V_post*(inv_V_prior*alpha_prior + kron(Z,invSIGMA)*Y(:));
    
    % check for stability of the VAR coefficients
    check=-1;
    while check<0        
        alpha = alpha_post + chol(V_post)'*randn(K*(const+K*p),1); % Draw of alpha
        A = reshape(alpha,K,const+K*p);                            % Draw of A
        Acomp = [A(:,2:end); eye(K*(p-1)) zeros(K*(p-1),K)];       % Companion matrix
        if (max(abs(eig(Acomp)))>1)==0                             % Check Eigenvalue of Companion matrix
            check = 1;
        end
    end
    
    % Posterior of SIGMA|alpha,Y ~ invWishard(inv(S_post),v_post)
    v_post = T + v_prior;
    resid_irep = Y - A*Z;
    S_post = S_prior + resid_irep*resid_irep';
    SIGMA = inv(wishrnd(inv(S_post),v_post));   % Draw of SIGMA

    % Store results if burn-in phase is passed
    if irep > nburn               
        A_draws(irep-nburn,:,:) = A;         % Store A
        SIGMA_draws(irep-nburn,:,:) = SIGMA; % Store SIGMA
    end
       
end
close(waitb);
toc; % Stop and display timer

%% INFERENCE ON POSTERIOR DRAWS
A_mean     = squeeze(mean(A_draws,1))     % Posterior mean of A
SIGMA_mean = squeeze(mean(SIGMA_draws,1)) % Posterior mean of SIGMA
A_OLS                                     % OLS estimate of A
SIGMA_OLS                                 % OLS estimate of SIGMA