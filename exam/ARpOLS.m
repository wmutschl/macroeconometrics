function OLS = ARpOLS(y,p,const,numAutoCorr,alph)
% =======================================================================
% OLS regression of AR(p) model:
% y_t = const + theta_1*y_{t-1} + ... + theta_p*y_{t-p} + u_t
% with u_t ~ WN(0,sigma_u)
% =======================================================================
% OLS = OLSmodel(y,p,const)
% -----------------------------------------------------------------------
% INPUT
%	- y     : dependent variable vector (T x 1)
%   - p     : number of lages
%   - const : 0 no constant; 1 constant; 2 constant and linear trend
% -----------------------------------------------------------------------
% OUTPUT
%	- OLS: structure including estimation results
%         - T_eff           : effective sample size used in estimation
%         - thetahat        : estimate of coefficients
%         - sig_uhat        : estimate of sigma_u
%         - sig_thetahat    : estimate of standard error of coefficients
%         - tstat           : t statistics, with alpha = 0.5
%         - pvalues         : p values of H_0: thetahat = 0
%         - theta_ci        : 95% confidence intervall for theta
% =======================================================================
% Willi Mutschler, October 2017
% willi@mutschler.eu
if nargin < 4
    numAutoCorr=10;
end
if nargin < 5
    alph = 0.05;
end

[T, K] = size(y);
T_eff = T-p; % effective sample size
Y = lagmatrix(y,1:p); % create matrix with lagged variables
if const==1 %constant
    Y = [ones(T_eff,1) Y((p+1):end,:)];
elseif const==2 % time trend and constant
    Y = [ones(T_eff,1) 1:T_eff' Y((p+1):end,:)];
else
    Y = Y((p+1):end,:);
end
y = y(p+1:end); % get rid of first p observations

YtYinv = inv(Y'*Y);
thetahat = YtYinv*(Y'*y);
yhat = Y*thetahat;
uhat = y - yhat;
utu = uhat'*uhat;

var_uhat = utu/(T_eff-p*K-const);
sig_uhat = sqrt(var_uhat);
var_thetahat = var_uhat*(diag(YtYinv));
sig_thetahat = sqrt(var_thetahat);

tstat = thetahat./sig_thetahat;
tcrit=-tinv(alph/2,T_eff-p*K-const);
pvalues = tpdf(tstat,T_eff-p*K-const);
theta_ci=[thetahat-tcrit.*sig_thetahat, thetahat+tcrit.*sig_thetahat];

r_k = nan(1,numAutoCorr);  % Initialize output vector
c0 = 1/(size(uhat,1))*(uhat' *uhat); % Compute variance
for h=1:numAutoCorr
    r_k(1,h) = 1/((size(uhat,1)-h)*c0) * (uhat(1+h:size(uhat,1),:)' * uhat(1:size(uhat,1)-h,:));
end
Q = size(uhat,1)*sum(r_k.*r_k);
Qp = chi2pdf(Q,h-p);

OLS.T_eff = T_eff;
OLS.thetahat = thetahat;
OLS.sig_uhat = sig_uhat;
OLS.sig_thetahat = sig_thetahat;
OLS.tstat = tstat;
OLS.pvalues = pvalues;
OLS.theta_ci = theta_ci;
OLS.uhat = uhat;
OLS.Q = Q;
OLS.Qp=Qp;
end