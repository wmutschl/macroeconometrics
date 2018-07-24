function [nlag,VAR] = LagOrderSelection(ENDO,EXOG,opt)
% =======================================================================
% Perform and display lag order selection tests
% =======================================================================
% LagOrderSelection(ENDO,EXOG,opt)
% -----------------------------------------------------------------------
% INPUTS
%   - ENDO : matrix [periods x number of endogenous variables]
%   - EXOG : matrix [periods x number of exogenous variables]
%   - opt  : structure with options, see load options section below
% -----------------------------------------------------------------------
% OUTPUTS
%   - nlag: number of lags choosen by the user, scalar
%   - VAR:  structure with VAR(nlag) estimation results
% -----------------------------------------------------------------------
% CALL
%   - EstimReducedForm.m: Perform Reduced Form VAR estimation with OLS
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

% Load options
const = opt.const;      % deterministic terms (scalar), values: 0 no constant, 1 constant, 2 constant+linear, 3 constant+linerar+quadratic
maxLag = opt.maxlags;   % number of maximum lags to consider (scalar)
nlag_ex = opt.nlag_ex;  % number of lags for exogenous variables (scalar)

if isempty(EXOG)==0 && (nlag_ex > maxLag && size(EXOG,1) > 0)
    maxLag = nlag_ex;    
end

% Initialize, +1 corresponds to no lags
AIC=nan(maxLag+1,1);
FPE=nan(maxLag+1,1);
BIC=nan(maxLag+1,1);
HQC=nan(maxLag+1,1);

% use same number of observations for all estimations
Y=ENDO(1+maxLag:end,:);
if isempty(EXOG)==0
    X = EXOG(1+maxLag:end,:);
end
[T,K] = size(Y); % T number of periods, K number of endogenous variables

opt.model = 1; % Estimate VAR model, not VECM
opt.coint_r = nan; % Cointegration rank is irrelevant here
for i=0:maxLag
    opt.nlag = i; % set lag
    p = sprintf('p%d',i);
    VAR.(p) = EstimReducedForm(Y,EXOG,opt); % Estimate reduced-form
    
    SIGU_ML = 1/T*VAR.(p).residuals'*VAR.(p).residuals;           % Maximum Likelihood estimation of covariance Matrix of reduced-form residuals
    AIC(i+1,:) = log(det(SIGU_ML))+2/T*(K*(K*i+const));           % Akaike information criteria
    BIC(i+1,:) = log(det(SIGU_ML))+log(T)/T*K*(K*i+const);        % Schwartz information criteria
    HQC(i+1,:) = log(det(SIGU_ML))+2*log(log(T))/T*K*(K*i+const); % Hannan-Quinn information criteria
end
% Store results and find minimal value of criteria
results = [transpose(0:maxLag) AIC BIC HQC];
pAIC = find(AIC == min(AIC))-1;
pBIC = find(BIC == min(BIC))-1;
pHQC = find(HQC == min(HQC))-1;

% Display summary of results
fprintf('******************************************************\n');
fprintf('*** OPTIMAL ENDOGENOUS LAGS FROM INFORMATION CRITERIA:\n');
fprintf('******************************************************\n');
disp(array2table(results,'VariableNames',{'Lag','AIC','BIC','HQC'}))
fprintf('  Optimal number of lags (searched up to %d lags of levels):\n',maxLag);
fprintf('  Akaike Info Criterion:    %d\n',pAIC);
fprintf('  Schwarz Criterion:        %d\n',pBIC);
fprintf('  Hannan-Quinn Criterion:   %d\n',pHQC);
fprintf('\n');

% Let the user decide on lag order and save corresponding VAR estimation
nlag = input('  Please set the lag order for the subsequent analysis: ');
VAR = VAR.(sprintf('p%d',nlag));

% DiagnosticTests(ENDO,EXOG,opt); % Diagnostic Tests against autocorrelation, nonnormality and ARCH effects in VAR residuals not yet
