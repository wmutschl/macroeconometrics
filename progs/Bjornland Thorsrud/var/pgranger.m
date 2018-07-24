function [ftest,fprob]=pgranger(y,x,constant,nlag,nvary,resid0)
% PURPOSE: Do Granger Causality test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% y = matrix (t x n), with dependent variables in VAR(nlag)
%
% x = matrix (t x m), with independent variables in VAR(nlag), where m is
% n*nlag + 1 (if constant). 
%
% constant = logical 
%
% nlag = number of lags
%
% nvary = n
%
% residu0 = matrix (t x n) form VAR(nlag) unrestricted regressions
%
% Output:
%
% ftest = matrix (n x n) with F-statistics 
%
% fprob = matrix (n x n) with p-values. A p-value less than e.g. 0.05
% indicates that we can reject the NULL of no-granger causality at the 5
% percent significance level. 
%
% Usage: 
%
% [ftest,fprob]=pgranger(y,x,constant,nlag,nvary,resid0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ftest=nan(nvary,nvary);                                        
fprob=nan(nvary,nvary);                                        

if constant
    xmat=x(:,2:end);
else
    xmat=x;
end;

nvarx=size(xmat,2);
nobs=size(y,1);

for j=1:nvary
    sigu=resid0(:,j)'*resid0(:,j);
    for r=1:nvary
        xtmp = [];
        for s=1:nvary;
            if s~=r
                xlag=xmat(:,s:nvary:nvarx);
                xtmp=[xtmp xlag];
            end;
        end;
        % we have an xtmp matrix that excludes 1 variable
        % add deterministic variables (if any) and constant term
        if constant
            xtmp=[ones(nobs,1) xtmp];               
        end;                
        
        % get ols residual vector
        b=xtmp\y(:,j); 
        etmp=y(:,j)-xtmp*b;
        sigr=etmp'*etmp;               
        % joint F-test for variables r
        ftest(j,r)=((sigr - sigu)/nlag)/(sigu/(nobs-nvarx)); 
        fprob(j,r)=1-fcdf(ftest(j,r),nlag,nobs-nvarx);                   
    end;
end;
