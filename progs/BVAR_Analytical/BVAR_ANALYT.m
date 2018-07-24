%--------------------------------------------------------------------------
% Bayesian estimation, prediction and impulse response analysis in VAR
% models. Dependent on your choice of forecasting, the VAR model is:
%
% Iterated forecasts:
%     Y(t) = A0 + Y(t-1) x A1 + ... + Y(t-p) x Ap + e(t)
%
% so that in this case there are p lags of Y (from 1 to p).
%
% Direct h-step ahead foreacsts:
%     Y(t+h) = A0 + Y(t) x A1 + ... + Y(t-p+1) x Ap + e(t+h)
%
% so that in this case there are also p lags of Y (from 0 to p-1).
%
% In any of the two cases, the model is written as:
%
%                   Y(t) = X(t) x A + e(t)
%
% where e(t) ~ N(0,SIGMA), and A summarizes all parameters. Note that we
% also use the vector a which is defined as a=vec(A). 
%
%--------------------------------------------------------------------------

tic;

%------------------------------LOAD DATA-----------------------------------
% Load Quarterly US data on inflation, unemployment and interest rate, 
% 1953:Q1 - 2006:Q3
load Yraw.dat;


% In any case, name the data you load 'Yraw', in order to avoid changing the
% rest of the code. Note that 'Yraw' is a matrix with T rows by M columns,
% where T is the number of time series observations (usually months or
% quarters), while M is the number of VAR dependent macro variables.
%----------------------------PRELIMINARIES---------------------------------
% Define specification of the VAR model
constant = 1;        % 1: if you desire intercepts, 0: otherwise 
p = 2;               % Number of lags on dependent variables
forecasting = 1;     % 1: Compute h-step ahead predictions, 0: no prediction
forecast_method = 0; % 0: Direct forecasts 
                     % 1: Iterated forecasts
h = 4;               % Number of forecast periods

% Set prior for BVAR model:
prior = 3;  % prior = 1 --> Noninformative Prior
            % prior = 2 --> Minnesota Prior
            % prior = 3 --> Natural conjugate Prior
  
%--------------------------DATA HANDLING-----------------------------------
% Get initial dimensions of dependent variable
[Traw M] = size(Yraw);

% The model specification is different when implementing direct forecasts,
% compared to the specification when computing iterated forecasts.
if forecasting==1
    if h<=0
        error('You have set forecasting, but the forecast horizon choice is wrong')
    end    

    % Now create VAR specification according to forecast method
    if forecast_method==0       % Direct forecasts
        Y1 = Yraw(h+1:end,:);
        Y2 = Yraw(2:end-h,:);
        Traw = Traw - h - 1;
    elseif forecast_method==1   % Iterated forecasts
        Y1 = Yraw;
        Y2 = Yraw;
    else
        error('Wrong choice of forecast_method')
    end
else
   Y1 = Yraw;
   Y2 = Yraw;
end
        
% Generate lagged Y matrix. This will be part of the X matrix
Ylag = mlag2(Y2,p); % Y is [T x M]. ylag is [T x (Mp)]

% Now define matrix X which has all the R.H.S. variables (constant, lags of
% the dependent variable and exogenous regressors/dummies)
if constant
    X1 = [ones(Traw-p,1) Ylag(p+1:Traw,:)];
else
    X1 = Ylag(p+1:Traw,:);  %#ok<UNRCH>
end

% Get size of final matrix X
[Traw3 K] = size(X1);

% Create the block diagonal matrix Z
Z1 = kron(eye(M),X1);

% Form Y matrix accordingly
% Delete first "LAGS" rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:); % This is the final Y matrix used for the VAR

% Traw was the dimesnion of the initial data. T is the number of actual 
% time series observations of Y and X
T = Traw - p;

%========= FORECASTING SET-UP:
% Now keep also the last "h" or 1 observations to evaluate (pseudo-)forecasts
if forecasting==1
    if forecast_method==0  % Direct forecasts, we only need to keep the 
        Y = Y1(1:end-1,:);                             % last observation
        X = X1(1:end-1,:);
        Z = kron(eye(M),X);
        T = T - 1;
    else              % Iterated forecasts, we keep the last h observations
        Y = Y1(1:end-h,:);
        X = X1(1:end-h,:);
        Z = kron(eye(M),X);
        T = T - h;
    end
else
    Y = Y1;
    X = X1;
    Z = Z1;
end

%--------------------------------PRIORS------------------------------------
% First get Ordinary Least Squares (OLS) estimators
A_OLS = inv(X'*X)*(X'*Y); % This is the matrix of regression coefficients
a_OLS = A_OLS(:);         % This is the vector of coefficients, i.e. it holds
                          % that a_OLS = vec(A_OLS)
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
SIGMA_OLS = SSE./(T-K);

%-----------------Prior hyperparameters for bvar model
% Define hyperparameters
if prior == 1 % Noninformtive
    % I guess there is nothing to specify in this case!
    % Posteriors depend on OLS quantities
elseif prior == 2 % Minnesota
    A_prior = 0*ones(K,M);   
    a_prior = A_prior(:);
    
    % Hyperparameters on the Minnesota variance of alpha
    a_bar_1 = 0.5;
    a_bar_2 = 0.5;
    a_bar_3 = 10^2;
    
    % Now get residual variances of univariate p-lag autoregressions. Here
    % we just run the AR(p) model on each equation, ignoring the constant
    % and exogenous variables (if they have been specified for the original
    % VAR model)
    sigma_sq = zeros(M,1); % vector to store residual variances
    for i = 1:M
        % Create lags of dependent variable in i-th equation
        Ylag_i = mlag2(Yraw(:,i),p);
        Ylag_i = Ylag_i(p+1:Traw,:);
        % Dependent variable in i-th equation
        Y_i = Yraw(p+1:Traw,i);
        % OLS estimates of i-th equation
        alpha_i = inv(Ylag_i'*Ylag_i)*(Ylag_i'*Y_i);
        sigma_sq(i,1) = (1./(T-p+1))*(Y_i - Ylag_i*alpha_i)'*(Y_i - Ylag_i*alpha_i);
    end
    
    % Now define prior hyperparameters.
    % Create an array of dimensions K x M, which will contain the K diagonal
    % elements of the covariance matrix, in each of the M equations.
    V_i = zeros(K,M);
    
    % index in each equation which are the own lags
    ind = zeros(M,p);
    for i=1:M
        ind(i,:) = constant+i:M:K;
    end
    for i = 1:M  % for each i-th equation
        for j = 1:K   % for each j-th RHS variable
            if constant==1
                if j==1
                    V_i(j,i) = a_bar_3*sigma_sq(i,1); % variance on constant                
                elseif find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(p^2); % variance on own lags           
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;                   
                        end
                    end
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));           
                end
            else
                if find(j==ind(i,:))>0
                    V_i(j,i) = a_bar_1./(p^2); % variance on own lags
                else
                    for kj=1:M
                        if find(j==ind(kj,:))>0
                            ll = kj;
                        end                        
                    end
                    V_i(j,i) = (a_bar_2*sigma_sq(i,1))./((p^2)*sigma_sq(ll,1));            
                end
            end
        end
    end
    
    % Now V is a diagonal matrix with diagonal elements the V_i
    V_prior = diag(V_i(:));  % this is the prior variance of the vector a  
    
    % SIGMA is equal to the OLS quantity
    SIGMA = SIGMA_OLS;
    
elseif prior == 3 % Normal-Wishart (nat conj)
    % Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
    A_prior = 0*ones(K,M);   
    a_prior = A_prior(:);
    V_prior = 10*eye(K);
    % Hyperparameters on inv(SIGMA) ~ W(v_prior,inv(S_prior))
    v_prior = M;
    S_prior = eye(M);
    inv_S_prior = inv(S_prior);
end
    
%============================ POSTERIORS ==================================
%==========================================================================
    
%--------- Posterior hyperparameters of ALPHA and SIGMA with Diffuse Prior
if prior == 1
    % Posterior of alpha|Data ~ Multi-T(kron(SSE,inv(X'X)),alpha_OLS,T-K)
    V_post = inv(X'*X);
    a_post = a_OLS;
    A_post = reshape(a_post,K,M);
    
    % posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
    S_post = SSE;
    v_post = T-K;
    
    % Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean = a_post;
    alpha_var = (1/(v_post - M - 1))*kron(S_post,V_post);
    
%--------- Posterior hyperparameters of ALPHA and SIGMA with Minnesota Prior
elseif prior == 2
    % ******Get all the required quantities for the posteriors       
    V_post = inv( inv(V_prior) + kron(inv(SIGMA),X'*X) );
    a_post = V_post * ( inv(V_prior)*a_prior + kron(inv(SIGMA),X'*X)*a_OLS );
    A_post = reshape(a_post,K,M);
     
    % In this case, the mean is a_post and the variance is V_post
   alpha_mean=a_post;
%--------- Posterior hyperparameters of ALPHA and SIGMA with Normal-Wishart Prior
elseif prior == 3
    % ******Get all the required quantities for the posteriors       
    % For alpha
    V_post = inv( inv(V_prior) + X'*X );
    A_post = V_post * ( inv(V_prior)*A_prior + X'*X*A_OLS );
    a_post = A_post(:);
    
    % For SIGMA
    S_post = SSE + S_prior + A_OLS'*X'*X*A_OLS + A_prior'*inv(V_prior)*A_prior - A_post'*( inv(V_prior) + X'*X )*A_post;
    v_post = T + v_prior;
    
    % Now get the mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean = a_post;
    alpha_var = (1/(v_post - M - 1))*kron(S_post,V_post); 
end

%======================= PREDICTIVE INFERENCE =============================
%==========================================================================

X_tplus1 = [1 Y(T,:) X(T,2:M*(p-1)+1)];

% As a point forecast use predictive mean 

Pred_mean = X_tplus1*A_post;


% Print some results
disp('The mean of alpha is in the vector alpha_mean')
disp('Its variance is in the vector alpha_var')
disp('Point forecast is Pred_mean')


toc;

