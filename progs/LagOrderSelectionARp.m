function nlag = LagOrderSelectionARp(ENDO,const,pmax)
% =======================================================================
% Perform and display lag order selection tests for AR(p) model, i.e.
% Akaike, Schwartz and Hannan-Quinn information criteria
% =======================================================================
% LagOrderSelection(ENDO,opt)
% -----------------------------------------------------------------------
% INPUTS
%   - ENDO  : data matrix [periods x 1]
%   - const : 0 no constant, 1 constant, 2 constant+linear trend. [scalar]
%   - pmax  : number of maximum lags to consider. [scalar]
% -----------------------------------------------------------------------
% OUTPUTS
%   - nlag: number of lags choosen by the user, scalar
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu

T = size(ENDO,1);     % Sample size
T_eff = T-pmax;       % Effective sample size used for all estimations, i.e.
                      % number of presample values set aside for estimation
                      % is determined by the maximum order pmax
% Initialize
AIC=nan(pmax,1);
SIC=nan(pmax,1);
HQC=nan(pmax,1);

% Construct regressor matrix and dependent variable, where the number of 
% presample values set aside for estimation is determined by pmax
if const==1     % constant
    YMAX = ones(T_eff,1);
elseif const==2 % constant and time trend
    YMAX = [ones(T_eff,1) transpose((pmax+1):T) ];
else
    YMAX = [];
end
YY = lagmatrix(ENDO,1:pmax); 
YMAX = [YMAX YY(pmax+1:T,:)];
y=ENDO(pmax+1:T,:);   

for p=1:pmax
    Y = YMAX(:,1:p+const);        % Data used in estimation
	thetahat = inv(Y'*Y)*Y'*y;    % OLS and ML estimator
    uhat = y - Y*thetahat;        % Residuals
    sigma2u=uhat'*uhat/T_eff;     % ML estimate of variance of errors
    np=p+const;                   % Number of freely estimated parameters
    logsigu2 = log(det(sigma2u));    
    AIC(p,:) = logsigu2 + 2/T_eff*np;                 % Akaike
    SIC(p,:) = logsigu2 + log(T_eff)/T_eff*np;        % Schwartz
    HQC(p,:) = logsigu2 + 2*log(log(T_eff))/T_eff*np; % Hannan-Quinn
end
% Store results and find minimal value of criteria
results = [transpose(1:pmax) AIC SIC HQC];
pAIC = find(AIC == min(AIC));
pBIC = find(SIC == min(SIC));
pHQC = find(HQC == min(HQC));

% Display summary of results
fprintf('******************************************************\n');
fprintf('*** OPTIMAL ENDOGENOUS LAGS FROM INFORMATION CRITERIA:\n');
fprintf('******************************************************\n');
disp(array2table(results,'VariableNames',{'Lag','AIC','BIC','HQC'}))
fprintf('  Optimal number of lags (searched up to %d lags):\n',pmax);
fprintf('  Akaike Info Criterion:    %d\n',pAIC);
fprintf('  Schwarz Criterion:        %d\n',pBIC);
fprintf('  Hannan-Quinn Criterion:   %d\n',pHQC);
fprintf('\n');

% Let the user decide on lag order
nlag = input('  Please set the lag order for the subsequent analysis: ');