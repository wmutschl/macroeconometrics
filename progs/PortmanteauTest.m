% =======================================================================
% Comparison of Portmanteau Tests For Residual Autocorrelation on inflation
% series for (i) AR(phat) where phat is determined by the AIC criteria and
% (ii) AR(1) model
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================
clearvars; clc; close all;
ENDO = xlsread('../data/gnpdeflator.xlsx',1); % load data
% computation of  inflation series
price_index = ENDO(:,3);
infl = log(price_index(2:end,:))- log(price_index(1:(end-1),:));

pmax = 12;                                    % set maximum number of lags
const = 0;                                    % model with constant
alph = 0.05;                                  % significance level
phat = LagOrderSelectionARp(infl,const,pmax); % compute criteria

OLSAR_phat = ARpOLS(infl,phat,const,alph);    % estimate AR(phat)
OLSAR_1    = ARpOLS(infl,1,const,alph);       % estimate AR(1)

% Compute portmanteu statistic
h = phat+10;          % maximum number of lags
u_p = OLSAR_phat.resid;
u_1 = OLSAR_1.resid;
T_p = size(u_p,1);
T_1 = size(u_1,1);
% initialize output vectors
rho_p = nan(1,h);
rho_1 = nan(1,h);
% compute variances
gam_p = 1/T_p*(u_p' *u_p);
gam_1 = 1/T_1*(u_1' *u_1);
% compute autocorrelations
for j=1:h
    rho_p(1,j) = 1/((T_p-j)*gam_p)*(u_p(1+j:T_p,:)'*u_p(1:T_p-j,:));
    rho_1(1,j) = 1/((T_1-j)*gam_1)*(u_1(1+j:T_1,:)'*u_1(1:T_1-j,:));
end
% compute test statistic
Q_p = T_p*sum(rho_p.*rho_p);
Q_1 = T_1*sum(rho_1.*rho_1);
% compute p values from chi2 distribution
Qpval_p = chi2pdf(Q_p,j-phat);
Qpval_1 = chi2pdf(Q_1,j-1);
% compare p values for AR(phat) and AR(1) model
disp([Qpval_p Qpval_1])