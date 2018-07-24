%%% Example 1: Monetary VAR
%%% 3 variables: real GDP, CPI, federal funds rate
%%% Sample period: 1970Q1 to 2011Q1
%%% Identification of 3 shocks with sign restrictions: 
%%% monetary policy shock, aggregate demand and aggregate supply shock
%1985-2006 data in first differences
% Algorithm follows Andrew Binning 2013, modified by
% Willi Mutschler, willi@mutschler.eu, December 2017

%% Set up model, restrictions and options
clearvars; clc;close all;
% Data Handling
Yraw = xlsread('../data/data.xlsx','MonPolData');
vnames   = {'GDP level','CPI','Federal Funds Rate'};
epsnames = {'Aggregate Demand','Aggregate Supply','Monetary Policy'};
IRFcumsum = [1 1 0];
nlag = 4;                 % Number of lags. [Scalar]
p = nlag;
const = 1;                % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend
StartValueMethod = 2;     %0: Use identity matrix, 1: use square root, 2: use cholesky as starting value
nsteps=21;                %horizon for impulse responses                                             
ndraws=200;               %number of solutions that match the sign restrictions
bmed=0.5;                 %impulse response percentile (median)
bupp=0.84;                %upper error band (percentile)
blow=0.16;                %lower error band (percentile)
it_print = 10;            % Print on the screen every "it_print"-th iteration

% Define Minnesota prior for BVAR model
hyperparams(1) = 0.5;   %tightness parameter for Minnesota prior on lags of own variable
hyperparams(2) = 0.5;   %tightness parameter for Minnesota prior on lags of other variable
hyperparams(3) = 100;     %tightness of prior on exogenous variables, i.e. constant term, trends, etc
% Gibbs-related preliminaries
nsave = 100;          % Final number of draws to save
nburn = 30000;          % Draws to discard (burn-in)
ntot = nsave+nburn;     % Total number of draws
it_print_bayes = 1000;        % Print on the screen every "it_print"-th iteration

% Sign restrictions on impact
Rsign  = [+1  +1   -1;
          +1  -1   -1;
          +1  nan  +1;
         ];

%% Estimate reduced-form model
%--------------------------DATA HANDLING-----------------------------------
[Traw,K] = size(Yraw);      % Get initial dimensions of dependent variable
Ylag = lagmatrix(Yraw,1:nlag); % Generate lagged Y matrix. This will be part of the Z matrix
% Now define matrix Z which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies), Also get rid of NA observations
if const == 0
    Z = transpose(Ylag(p+1:Traw,:));  
elseif const == 1
    Z = transpose([ones(Traw-p,1) Ylag(p+1:Traw,:)]);
elseif const == 2
    Z = transpose([ones(Traw-p,1) transpose((p+1):Traw) Ylag(p+1:Traw,:)]);
end
Y = transpose(Yraw(p+1:Traw,:)); % Dependent variable in i-th equation, get rid of observations
[totcoeff,T] = size(Z); % Get size of final matrix Z
ZZt = Z*Z'; % auxiliary product

%-----------------------------INITIALIZATION--------------------------------
% Get OLS estimators
A_OLS = (Y*Z')/ZZt; 
a_OLS = A_OLS(:);
resid_OLS = Y - A_OLS*Z;
SSE = resid_OLS*resid_OLS';
SIGMA_OLS = SSE./(T-K*p-const);

% Initialize Bayesian posterior parameters using OLS values
alpha = a_OLS;     % This is the single draw from the posterior of alpha
A = A_OLS;     % This is the single draw from the posterior of ALPHA
SSE_Gibbs = SSE;   % This is the single draw from the posterior of SSE
SIGMA = SIGMA_OLS; % This is the single draw from the posterior of SIGMA

%--------------------------PRIOR SPECIFICATION-----------------------------------
% Get standard versions of Minnesota
[a_prior, V_prior, inv_V_prior, v_prior, S_prior, inv_S_prior] = BVARMinnesotaPrior(Yraw,const,nlag,hyperparams);

%========================== Start Sampling ================================
%==========================================================================
% Storage space for posterior draws
%A_draws = zeros(nsave,K,totcoeff);
%SIGMA_draws = zeros(nsave,K,K);
nvars =K;
I_K    = eye(nvars);
I_KK   = eye(nvars^2);
J = [I_K zeros(nvars,nvars*(nlag-1))];
%selSign = find(isnan(Rsign)==0); % Index for sign restrictions on impact

IRFvals = zeros(nvars,nsteps,nsave*ndraws,nvars); % Contains impulse resonse functions
% 1st dimension = variable
% 2nd dimension = time
% 3rd dimension = rotation draw
% 4th dimension = shock

for irep = 1:ntot  %Start the Gibbs "loop"
    %Posterior of ALPHA|SIGMA,Data ~ N(a_post,V_post)
    invSIGMA = inv(SIGMA);    
    V_post = (inv_V_prior+kron(ZZt,invSIGMA))\eye(size(inv_V_prior,1));
    a_post = V_post*(inv_V_prior*a_prior + kron(Z,invSIGMA)*Y(:));
    %check for stability of the VAR
    check=-1;
    while check<0        
        alpha = a_post + chol(V_post)'*randn(K*(const+K*p),1); % Draw of alpha
        A = reshape(alpha,K,const+K*p); % Draw of A
        Acomp = [A(:,2:end); eye(K*(p-1)) zeros(K*(p-1),K)];
        if (max(abs(eig(Acomp)))>1)==0
            check = 1;
        end
    end
    % Posterior of SIGMA|ALPHA,Data ~ iW(inv(S_post),v_post)
    v_post = T + v_prior;
    resid_irep = Y - A*Z;
    S_post = S_prior + resid_irep*resid_irep';
    SIGMA = inv(wishrnd(inv(S_post),v_post));% Draw SIGMA
    % Estimate Structural Model    
    if irep > nburn 
        disp(irep-nburn)
        B0inv_0 = chol(SIGMA,'lower'); % Use Cholesky decomposition of vcov of reduced form
        accepted_draws = 1;
        discarded_draws = 0;
        while accepted_draws < ndraws+1
            % Draw an orthogonal matrix Q from set set of all orthogonal matrices
            W = randn(nvars,nvars);
            [Q,R] = qr(W);
            % Normalization: make sure that diagonal of upper triangular R is positive
            for i = 1:nvars
                if R(i,i)<0 % if value on diagonal negative
                    Q(:,i) = -Q(:,i); % reverse all signs in the i-th column of Q
                end
            end
            B0inv = B0inv_0*Q; % New draw for B0inv    

            % Check sign restrictions
            for jj = 1:nvars
                iota = I_K(:,jj);
                chk = B0inv*iota;
                THETA = zeros(nvars*nlag,nsteps); % initialize impulse response function used in companion form
                THETA(1:nvars,1) = B0inv*iota;    % initial impact of shock jj

                sr_index = ~isnan(Rsign(:,jj));
                tmp = sign(chk(sr_index)) - Rsign(sr_index,jj);
                % if sign restrictions are not fullfilled, return to beginning of while loop
                if any(tmp~=0)                
                    jj = 0;
                    discarded_draws = discarded_draws + 1; % Count discarded draws
                    break 
                end
                % Else compute IRF for nsteps for companion form
                for ii = 2:nsteps
                    THETA(:,ii) = Acomp*THETA(:,ii-1);                
                end
                % Save IRF of shock j into output structure, we only need first K rows of companion form
                IRFvals(:,:,(irep-nburn-1)*nsave+accepted_draws,jj) = THETA(1:nvars,:);
            end

            if jj == nvars
                accepted_draws = accepted_draws + 1; % Count accepted draws
            end
        end
    end
       
end %end the main Gibbs for loop

%====================== End Sampling Posteriors ===========================

% Display median response of impulse response functions that satisfy restrictions
    figure('Name','Median Response And Percentiles');
    countplots = 1;
    for ishock = 1:nvars
        for ivar = 1:nvars
            subplot(nvars,nvars,countplots);
            IRFS = squeeze(IRFvals(ivar,1:nsteps,:,ishock));
            if IRFcumsum(ivar) == 1
                IRFS = cumsum(IRFS,1);
            end        
            IRFSmed = median(IRFS,2);
            IRFSup  = prctile(IRFS,bupp*100,2);
            IRFSlow  = prctile(IRFS,blow*100,2);
            plot(0:1:nsteps-1, IRFSmed ,'b','linewidth',2);
            hold on;
            plot(0:1:nsteps-1, [IRFSlow IRFSup] ,'--r');
            grid;
            title(vnames{ivar})
            ylabel([epsnames{ishock}, 'Shock'])
            countplots = countplots + 1;
        end
    end
