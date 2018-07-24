function ReducedForm = EstimReducedForm(ENDO,EXOG,opt)
% =======================================================================
% Perform reduced-form estimation: VAR model is estimated via OLS, 
% VECM model is estimated via reduced-rank maximum likelihood a la Johansen
% =======================================================================
% ReducedForm = EstimReducedForm(ENDO,EXOG,opt)
% -----------------------------------------------------------------------
% INPUTS
%	- ENDO : Matrix of endogenous variables. [periods x number of endogenous variables]
%	- EXOG : Matrix of exogenous variables. [periods x number of exogenous variables]
%   - opt  : Structure with options, see load options section below
% -----------------------------------------------------------------------
% OUTPUTS
%   - ReducedForm: structure including reduced form estimation results with the following fields:
%       - ENDO       : Matrix of endogenous variables. [periods x nvars]
%       - EXOG       : Matrix of exogenous variables. [periods x nvars_ex]
%       - AMATt      : Matrix of coefficients of VAR representation. [nvars x nvars*nlags]
%       - Acomp      : Companion matrix of VAR representation. [nvars*nlags x nvars*nlags]
%       - maxEig     : Maximum Eigenvalue of companion matrix. scalar
%       - residuals  : Reduced-form residuals. [periods x 1]
%       - SIGMAUHAT  : Estimated covariance matrix of reduced-form residuals. [nvars x nvars]
%       - Y          : Effective data matrix of endogenous variables used for OLS estimation
%       - X          : Effective data matrix of lagged endogenous variables used for OLS estimation
%       - X_EX       : Matrix of all exogenous variables used for OLS estimation
%       - nobs       : Sample size used in estimation
%       - paramVals  : Structure of parameter estimates with field names corresponding to the parameter names in paramNames (only VECM)
%       - paramNames : Cell vector of parameter names (only VECM)
%       - coint_r    : Specified cointegration rank, scalar (only VECM)

% -----------------------------------------------------------------------
% CALLS
%   - VARmakexy.m to create independent vector and lagged dependent matrix
%   - VARmakelags.m to build a matrix with lagged values of data
%   - jcitest.m to estimate VECM with reduced-rank maximum likelihood
%   - vec2var.m to convert VECM model to VAR model
% =======================================================================
% This function is inspired by Ambrogio Cesa Bianchi's VARmodel.m and 
% James P. LeSage's vare.m. VECM estimation is added.
% -----------------------------------------------------------------------
% Willi Mutschler, March 2017
% willi@mutschler.eu


%% Load options
%===============================================
[nobs, nvars] = size(ENDO); % nobs number of observations, nvars number of endogenous variables
nlag = opt.nlag;            % number of lags for endogenous variables, i.e. VAR(nlag) or VECM(nlag-1) is estimated
nlag_ex = opt.nlag_ex;      % number of lags for exogenous variables
const = opt.const;          % 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend, and trend^2
JOHtest = opt.JOHtest;         % String indicating the type of cointegratoin test to be performed in jcitest.m. Values are 'trace' or 'maxeig'
JOHmodel = opt.JOHmodel;       % String specifying the form of the deterministic components of the VEC(nlag-1) model. Values are 'H2','H1*','H1','H*','H', see jcitest.m for details
if isfield(opt,'JOHCointTest')
    JOHCointTest = opt.JOHCointTest; % Display cointegration test results
else
    JOHCointTest = 0; % Only estimate VECM model
end
% if isfield(opt,'coint_r')
%     coint_r = opt.coint_r; % 
% end

%% Create data matrices
%===============================================
% Check if there are exogenous variables 
if isempty(EXOG)==0
    [nobs2, nvar_ex] = size(EXOG); %nobs2 number of periods of exogenous variables, nvar_ex number of exogenous variables
    % Check that ENDO and EXOG are conformable
    if (nobs2 ~= nobs)
        error('var: nobs in EXOG-matrix not the same as y-matrix');
    end
    clear nobs2
    % Check if there is lag order of EXOG, otherwise set it to 0
    if isempty(nlag_ex)
        nlag_ex = 0;
    end
    if opt.model == 2 % Note that for VECM models jcitest does not support other exogenous variables
        error('jcitest does not support exogenous variables, discarding them.');
    end
else
    nvar_ex = 0;
    nlag_ex = 0;
    EXOG = [];
end

nobse     = nobs - max(nlag,nlag_ex);   % number of periods used in estimation
ncoeff    = nvars*nlag;                 % number of endogenous coefficients estimated in VAR model
ncoeff_ex = nvar_ex*(nlag_ex+1);        % number of exogenous coefficients estimated
ntotcoeff = ncoeff + ncoeff_ex + const; % total number of coefficients to be estimated in VAR
    
% Create independent vector and lagged dependent matrix
[Y, X] = VARmakexy(ENDO,nlag,const);

% Create (lagged) exogenous matrix
if nvar_ex>0
    X_EX  = VARmakelags(EXOG,nlag_ex);
    if nlag == nlag_ex
        X = [X X_EX];
    elseif nlag > nlag_ex
        diff = nlag - nlag_ex;
        X_EX = X_EX(diff+1:end,:);
        X = [X X_EX];
    elseif nlag < nlag_ex
        diff = nlag_ex - nlag;
        Y = Y(diff+1:end,:);
        X = [X(diff+1:end,:) X_EX];
    end
else
    X_EX = [];
end

    
%% Compute the matrix of coefficients & Variance Covariance Matrix
%===============================================
if opt.model == 1 %VAR model
    switch opt.EstMethodVAR
        case 'OLS'
            AMATt = (X'*X)\(X'*Y);
            UHAT = Y-X*AMATt;
            SIGMAUHAT = (1/(nobse-ntotcoeff))*(UHAT'*UHAT); % adjusted for number of estimated coefficients
        case 'ML'
            AMATt = (X'*X)\(X'*Y);
            UHAT = Y-X*AMATt;
            SIGMAUHAT = (1/(nobse))*(UHAT'*UHAT); % equivalent to maximum likelihood estimation
    end
    AMAT = AMATt'; 
    % Companion matrix
    if nlag >0
        Acomp = [AMAT(:,1+const:nvars*nlag+const); eye(nvars*(nlag-1)) zeros(nvars*(nlag-1),nvars)];
    else
        Acomp = [];
    end
elseif opt.model == 2 % VECM model
    if JOHCointTest == 1
        [~,~,~,~,JOH_mles] = jcitest(ENDO,'model',JOHmodel,'lags',nlag-1,'test',JOHtest);
        % Let the user specify the cointegration rank
        coint_r = input('Please set the cointegration rank for the subsequent analysis: ');
    else
        warning off;
        [~,~,~,~,JOH_mles] = jcitest(ENDO,'model',JOHmodel,'lags',nlag-1,'test',JOHtest,'Display','off');
        warning on;
        coint_r = opt.coint_r;
    end
    MLES = JOH_mles.(sprintf('r%d',coint_r)); % Save ML estimation results for VECM(nlag-1) and cointegration rank coint_r
    % Transform VECM to get corresponding VAR representation
    if nlag > 1
        % Values of parameter matrices
        for i = 1:(nlag-1)
            VECM_mat{i} = MLES.paramVals.(sprintf('B%d',i));
        end
        % Value of constant term depending on deterministic terms in VECM model
        VECM_C = MLES.paramVals.A*MLES.paramVals.B';
        switch opt.JOHmodel
            case 'H2'
                VECM_D = zeros(nvars,1);
            case 'H1*'
                VECM_D = MLES.paramVals.A*MLES.paramVals.c0;
            case 'H1'
                VECM_D = MLES.paramVals.A*MLES.paramVals.c0 + MLES.paramVals.c1;
            case 'H*'
                VECM_D = [MLES.paramVals.A*MLES.paramVals.c0 + MLES.paramVals.c1, MLES.paramVals.A*MLES.paramVals.d0];
            case 'H'
                VECM_D = [MLES.paramVals.A*MLES.paramVals.c0 + MLES.paramVals.c1, MLES.paramVals.A*MLES.paramVals.d0 + MLES.paramVals.d1];
        end
        AMAT = [VECM_D cell2mat(vec2var(VECM_mat,VECM_C))];
        % Companion matrix
        Acomp = [AMAT(:,1+const:nvars*nlag+const); eye(nvars*(nlag-1)) zeros(nvars*(nlag-1),nvars)];            
        AMATt =transpose(AMAT);
    else
        Acomp = [];            
    end
    % Store additional results for VECM model
    ReducedForm.paramVals  = MLES.paramVals;
    ReducedForm.paramNames = MLES.paramNames;
    ReducedForm.coint_r = coint_r; %cointegration rank
    UHAT  = MLES.res; %residuals
    SIGMAUHAT  = MLES.EstCov; % estimated maximum likelihood covariance matrix of residuals    
end

%% Save other results
%===============================================

ReducedForm.ENDO = ENDO;
ReducedForm.EXOG = EXOG;
ReducedForm.AMATt = AMATt;
ReducedForm.Acomp = Acomp;
ReducedForm.maxEig = max(abs(eig(Acomp)));
ReducedForm.residuals = UHAT;
ReducedForm.SIGMAUHAT = SIGMAUHAT;
ReducedForm.Y = Y;
ReducedForm.X = X;
ReducedForm.X_EX = X_EX;
ReducedForm.nobs = nobse;

