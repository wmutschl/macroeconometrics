function result=estimateVAR(y,x,nvary,nlag)

nobs=size(y,1);
nvarx=size(x,2);
% OLS
coeff=(x\y)';
beta=coeff(:,1:nvary*nlag);
alfa=coeff(:,nvary*nlag:end);
yhat=(coeff*x')';
resid=y-yhat;
% ML estimate of covaraince of residuals
Sigma2=(resid'*resid)/nobs;
% Compute some statistics
sige=sqrt(diag((resid'*resid)/(nobs-nvarx)));
bstd=nan(nvary,nvarx);            
sigma2hat=sige(:).^2;
invx=diag(inv(x'*x));
for i=1:nvary 
    bstd(i,:)=sqrt(sigma2hat(i)*invx);            
end;
R2=nan(nvary,1);
m0=eye(nobs)-1/nobs.*ones(nobs);   
ee=diag(resid'*resid);                        
for i=1:nvary                                                 
    R2(i)=1-ee(i)/(y(:,i)'*m0*y(:,i));                        
end;
% Steady state of VAR
% Could alternatively do this by the companion form
%[betac,alfac]=varGetCompForm(beta,alfa,nlag,nvary);
% inv(eye(nlag^2)-betac)*alfac ... and select the nvary first elements
As=0;
for j=1:nlag
    As=As+beta(:,(j-1)*nvary+1:j*nvary);
end
SS=(eye(nvary)-As)\alfa;          
% Check stationarity
betaCbb=varGetCompForm(beta,[],nlag,nvary);
if any(abs(real(eig(betaCbb)))>=1)
    stabilityCondition=1;
else
    stabilityCondition=0;
end;
%% Put into result struct
result.alfa=alfa;
result.beta=beta;
result.Sigma2=Sigma2;
result.bstd=bstd;
result.R2=R2;
result.SS=SS;
result.stabilityCondition=stabilityCondition;
result.nobs=nobs;
result.nvarx=nvarx;
result.nvary=nvary;
result.nlag=nlag;
result.yhat=yhat;
result.y=y;
result.x=x;