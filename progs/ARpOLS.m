function OLS = ARpOLS(y,p,const,alph)
% =======================================================================
% OLS regression of AR(p) model:
% y_t = c + d*t + theta_1*y_{t-1} + ... + theta_p*y_{t-p} + u_t
% with white noise u_t ~ (0,sigma_u)
% =======================================================================
% OLS = OLSmodel(y,p,const,alph)
% -----------------------------------------------------------------------
% INPUT
%   - y     : data vector of dimension T. [Tx1]
%   - p     : number of lags. [scalar]
%   - const : 1 constant; 2 constant and linear trend in model. [scalar]
%	- alpha : significance level for t statistic and p value. [scalar]
% -----------------------------------------------------------------------
% OUTPUT
%	- OLS: structure including estimation results
%     - T_eff        : effective sample size used in estimation. [scalar]
%     - thetahat     : estimate of coefficients. [(const+p)x1]
%     - sig_thetahat : estimate of standard error of coefficients. [(const+p)x1]
%     - tstat        : t statistics. [(const+p)x1]
%     - pvalues      : p values of H_0: thetahat = 0. [(const+p)x1]
%     - sig_uhat     : estimate of standard deviation of error term u. [scalar]
%     - theta_ci     : (1-alph)% confidence intervall for theta given 
%                      significance level alph. [(const+p)x2]
%     - resid        : residuals. [T_eff x 1]
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================

T = size(y,1);            % sample size
T_eff = T-p;              % effective sample size used in estimation
Y = lagmatrix(y,1:p);     % create matrix with lagged variables
if const==1               % add constant term
    Y = [ones(T,1) Y];
elseif const==2           % add constant term and time trend
    Y = [ones(T,1) transpose(1:T) Y];
end
Y = Y((p+1):end,:);       % get rid of initial p observations
y = y(p+1:end);           % get rid of initial p observations

YtYinv = inv(Y'*Y);
thetahat = YtYinv*(Y'*y); % OLS estimator of coefficients
yhat = Y*thetahat;        % predicted values
uhat = y - yhat;          % residuals
utu = uhat'*uhat;       

var_uhat = utu/(T_eff-p-const);         % variance of error term
sig_uhat = sqrt(var_uhat);              % standard deviaiont of error term
var_thetahat = var_uhat*(diag(YtYinv)); % variance of coefficients
sig_thetahat = sqrt(var_thetahat);      % standard error of coefficients

tstat = thetahat./sig_thetahat;         % t-statistics
tcrit=-tinv(alph/2,T_eff-p-const);      % critical value
pvalues = tpdf(tstat,T_eff-p-const);    % p-value
% confidence interval
theta_ci=[thetahat-tcrit.*sig_thetahat, thetahat+tcrit.*sig_thetahat];

% Store into output structure
OLS.T_eff        = T_eff;
OLS.thetahat     = thetahat;
OLS.sig_uhat     = sig_uhat;
OLS.sig_thetahat = sig_thetahat;
OLS.tstat        = tstat;
OLS.pvalues      = pvalues;
OLS.theta_ci     = theta_ci;
OLS.resid        = uhat;
end