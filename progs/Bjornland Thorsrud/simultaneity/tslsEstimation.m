function results=tslsEstimation(y,y1,x1,xall)
% PURPOSE: Do Two-stage least-squares regression (TSLS or 2SLS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model:
%
% y1 = c_1*xall + u
% y  = b_1*y1 + b_2*x1 + e; Equation of interest
%
% Input:
%
% y = vector (n x 1). Dependent variable
%
% y1 = matrix (n x g). Endogenous variables
%
% x1 = matrix (n x k). Exogenous variables for the dependent equation
%
% xall = matrix (n x (k+i)). Exogenous variables + all instruments not
% included in the equation for the dependent variable, i.e. [x1 instruments]
%
% Output: 
%
% result = structure, with the following fields:
% results.meth  = string: 'tsls';
% results.bhat  = bhat estimates, vector ((g+k) x 1)
% results.bstd  = bstd estimates, vector ((g+k) x 1)
% results.tstat = t-statistics, vector ((g+k) x 1) 
% results.yhat  = yhat predicted values, vector (n x 1)
% results.resid = residuals, vector (n x 1)
% results.sige  = e'*e/(n-k), scalar
% results.rsqr  = rsquared, scalar
% results.rbar  = rbar-squared, scalar
% results.dw    = Durbin-Watson Statistic, scalar
% results.nobs  = number of observation, n
% results.nendog = # of endogenous
% results.nexog  = # of exogenous
% results.nvar   = results.nendog + results.nexog
% results.y      = y data vector, vector (n x 1)
%
% Usage: 
%
% result=tsls(y,y1,x1,x2)
%
% Note: A constant should be included in both x1 and x2. The ordering of
% the paramter estimates follows the ordering of the variables in the input
% and as dispalyed under Model above. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a slightly modified version of the code written by:
%
% James P. LeSage, Dept of Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin~=4; 
    error([mfilename ':input: Wrong # of arguments to ' mfilename']); 
end;
results.meth='tsls';
[nobs1,g]=size(y1);
[nobs2,k]=size(x1);
[nobs3,l]=size(xall);
results.nendog=g; 
results.nexog=k; 
results.nvar=k+g;
if nobs1==nobs2;
    if nobs2==nobs3
        nobs=nobs1;
    end;
else
    error([mfilename ':input: # of observations in yendog, xexog, xall not the same']);     
end;
results.y=y;
results.nobs=nobs;
% xall contains all explanatory variables
% x1 contains exogenous
% y1 contains endogenous
xapxa=inv(xall'*xall);
% form xpx
xpx=[y1'*xall*xapxa*xall'*y1     y1'*x1
       x1'*y1                      x1'*x1];
xpy=[y1'*xall*xapxa*xall'*y
       x1'*y                  ];
xpxi=inv(xpx);                
results.beta=xpxi*xpy;             % bhat
results.yhat=[y1 x1]*results.beta; % yhat
results.resid=y-results.yhat;     % residuals
sigu=results.resid'*results.resid;
results.sige=sigu/(nobs-k-g);       % sige
results.bstd=sqrt(results.sige*(diag(xpxi)));
results.tstat=results.beta./results.bstd;
ym=y-ones(nobs,1)*mean(y);
rsqr1=sigu;
rsqr2=ym'*ym;
results.rsqr=1.0-rsqr1/rsqr2; % r-squared
rsqr1=rsqr1/(nobs-k-g);
rsqr2=rsqr2/(nobs-1.0);
results.rbar=1-(rsqr1/rsqr2); % rbar-squared
ediff=results.resid(2:nobs)-results.resid(1:nobs-1);
results.dw=(ediff'*ediff)/sigu; % durbin-watson