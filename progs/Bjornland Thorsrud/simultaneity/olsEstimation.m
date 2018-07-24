function results=olsEstimation(y,x)
% PURPOSE: Do OLS estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(y,1);

beta_ols=x\y;
% some statistics
yhat=x*beta_ols;
e_ols=y-yhat;
sige=(diag((e_ols'*e_ols)/(N-size(x,2))));
invx=diag(inv(x'*x));
bstd_ols=sqrt(sige*invx);            

ym=y-ones(N,1)*mean(y);
rsqr1=e_ols'*e_ols;
rsqr2=ym'*ym;
rsqr=1.0-rsqr1/rsqr2; % r-squared

% Put studd into results struct:
results.beta=beta_ols;
results.bstd=bstd_ols;
results.tstat=beta_ols./bstd_ols;
results.rsqr=rsqr;
results.sige=sige;