function result=chapter8_estimateVARgdpAndUnemp()
% PURPOSE: Do Cholesky and Blanchard & Quah identification of structural VAR
% in GDP and unemployment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The code is not written for computational efficiency or
% elegance. The book: 
%
% "Applied Time Series for Macroeconomics"
% Gyldendal Akademisk 2014
% by Hilde C. Bjørnland and Leif A. Thorsrud 
%
% provides details. Please refer to the book if the code(s) are used for 
% research of commercial purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings

% Number of lags in VAR
nlag=2;
% Number of bootstraps (to simulate confidence bands)
draws=1000;
% Impulse-response horizon
nirf=21;

%% Load data

[data,txt]=xlsread('./blanchardquaData.xlsx','data');

y1=data(:,2);
y2=data(:,3);
dates=data(:,1);

Y=[y1 y2];
[T,nvary]=size(Y);
ynames={'GDP','U'};

%% Estimate by OLS and compute some statistics using lag length from above

X=[ones(T,1) latMlag(Y,nlag)];
y=Y(1+nlag:end,:);
x=X(1+nlag:end,:);

nobs=size(y,1);
nvarx=size(x,2);

% OLS
coeff=(x\y)';
beta=coeff(:,2:end);
alfa=coeff(:,1);
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

% Add to results struct
result.SS=SS;

%% Do a residual bootstrap of VAR 
[betac,alfac]=varGetCompForm(beta,alfa,nlag,nvary);
% bootstrapped paramete values
[alfab,betab,Sigma2b]=varSimulation(alfac,betac,resid,nlag,draws);

%% Compute MSE, IRFs and VDC using Cholesky and long run identification

% Point estimate mse
mse=getMse(betac,Sigma2,nirf);

% get comp form of simulated parameters
[betacb,alfacb]=varGetCompForm(betab,alfab,nlag,nvary);

% Cholesky identification
irf_ch=varImpulseResponsFast(betac,nvary,nlag,nirf,chol(Sigma2,'lower'));
vdc_ch=getVdcFromIrf(irf_ch);
%...get for different draws (i.e., simulate)
[irfb_ch,vdcb_ch]=deal(nan(nvary,nvary,nirf,draws));
for d=1:draws
    irfb_ch(:,:,:,d)=varImpulseResponsFast(betacb(:,:,d),nvary,nlag,nirf,chol(Sigma2b(:,:,d),'lower'));
    vdcb_ch(:,:,:,d)=getVdcFromIrf(irfb_ch(:,:,:,d));
end;

% Long run identification
A0=longRunIdentification(betac,Sigma2);
irf_lo=varImpulseResponsFast(betac,nvary,nlag,nirf,A0);
vdc_lo=getVdcFromIrf(irf_lo);
%...get for different draws (i.e., simulate)
[irfb_lo,vdcb_lo]=deal(nan(nvary,nvary,nirf,draws));
for d=1:draws
    A0=longRunIdentification(betacb(:,:,d),Sigma2b(:,:,d));
    irfb_lo(:,:,:,d)=varImpulseResponsFast(betacb(:,:,d),nvary,nlag,nirf,A0);
    vdcb_lo(:,:,:,d)=getVdcFromIrf(irfb_lo(:,:,:,d));
end;

% Add to result struct()
result.mse=mse;
result.irf_ch=irf_ch;
result.vdc_ch=vdc_ch;
result.irf_lo=irf_lo;
result.vdc_lo=vdc_lo;


%% Plot impulse responses
alphaSign=0.1;

% make level of GDP growth...
irf_ch(1,1,:)=cumsum(irf_ch(1,1,:),3);
irf_ch(1,2,:)=cumsum(irf_ch(1,2,:),3);
irf_lo(1,1,:)=cumsum(irf_lo(1,1,:),3);
irf_lo(1,2,:)=cumsum(irf_lo(1,2,:),3);
% ...and for the draws
irfb_ch(1,1,:,:)=cumsum(irfb_ch(1,1,:,:),3);
irfb_ch(1,2,:,:)=cumsum(irfb_ch(1,2,:,:),3);
irfb_lo(1,1,:,:)=cumsum(irfb_lo(1,1,:,:),3);
irfb_lo(1,2,:,:)=cumsum(irfb_lo(1,2,:,:),3);

% Get quantiles (important, do this after level correction)
irf_ch_q=getQuantiles(irfb_ch,alphaSign);
irf_lo_q=getQuantiles(irfb_lo,alphaSign);

%% Make impulse response figures

% Cholesky id
makeIRFfigure(permute(irf_ch(1,1,:),[3 1 2]),permute(irf_ch_q(1,1,:,[1 3]),[3 4 1 2]),{'GDP','GDP'},'ch');
makeIRFfigure(permute(irf_ch(1,2,:),[3 1 2]),permute(irf_ch_q(1,2,:,[1 3]),[3 4 1 2]),{'GDP','U'},'ch');
makeIRFfigure(permute(irf_ch(2,1,:),[3 1 2]),permute(irf_ch_q(2,1,:,[1 3]),[3 4 1 2]),{'U','GDP'},'ch');
makeIRFfigure(permute(irf_ch(2,2,:),[3 1 2]),permute(irf_ch_q(2,2,:,[1 3]),[3 4 1 2]),{'U','U'},'ch');
% Long run id
makeIRFfigure(permute(irf_lo(1,1,:),[3 1 2]),permute(irf_lo_q(1,1,:,[1 3]),[3 4 1 2]),{'GDP','GDP'},'lo');
makeIRFfigure(permute(-irf_lo(1,2,:),[3 1 2]),permute(-irf_lo_q(1,2,:,[1 3]),[3 4 1 2]),{'GDP','U'},'lo'); % - such that we look at a positive supply shock
makeIRFfigure(permute(irf_lo(2,1,:),[3 1 2]),permute(irf_lo_q(2,1,:,[1 3]),[3 4 1 2]),{'U','GDP'},'lo');
makeIRFfigure(permute(-irf_lo(2,2,:),[3 1 2]),permute(-irf_lo_q(2,2,:,[1 3]),[3 4 1 2]),{'U','U'},'lo'); % - such that we look at a positive supply shock

