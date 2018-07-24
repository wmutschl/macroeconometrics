function chapter2_estimateAR()
%% PURPOSE: Estimate AR(p) models for GDP and inflation using OLS. Lags etc. 
% based on information criteria
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

%% 1) Settings

% Maximum number of lags to consider in model
P = 8;

% Minimum number of observations when computing estimates using an
% expanding estimation window (parameter stabiity) 
windowSize=80;

variables={'gdp','cpi'};

%% 2) Load data. Construct growth rate
for i=1:numel(variables)


    if i==1

        % time index settings
        eStart = 1981.01;   % first period for estimation 
        eEnd = 2013.01;    

        [data,~] = xlsread('./gdpNorway.xlsx','gdp');
        dates = data(:,1);
        yLevel = data(:,2);
        % take log and get the log qoq diff (scaled by 100 to interpret as
        % percent)
        yLevelLog = log(yLevel).*100;
        y = yLevelLog(5:end,1)-yLevelLog(1:end-4,1);        
        % correct dates for one lost observation due to diff
        dates = dates(5:end);

    else

        % time index settings
        eStart = 1996.01;   % first period for estimation 
        eEnd = 2013.02;    
               
        data=xlsread('./infNorway.xlsx','m')';
        % Data already on yoy form
        y=data(:,2);
        dates =data(:,1);

    end;


    %% 3) Do aic and bic test to choose lag length. 

    %Do this on whole sample 
    [aicTest,bicTest] = deal(nan(P,1)); % empty output variables

    stIdx = find(dates==eStart);        % start index for estiamtion
    enIdx = find(dates==eEnd);      % end index for estimation
    for p=1:P    
        [yt,xt] = constructArModel(y,p,stIdx,enIdx);
        [~,~,errors] = estimateArModel(yt,xt);    
        bicTest(p,1) = BIC(errors,p);
        aicTest(p,1) = AIC(errors,p);        
    end;
    % Get number of lags based on aicTest
    [~,li]=min(aicTest(:,1));
    nlag=li;

    %% Re-estimate model using this lag
    [yt,xt] = constructArModel(y,nlag,stIdx,enIdx);
    [yhat,beta,resid] = estimateArModel(yt,xt);    
    nobs=size(resid,1);
    nvarx=nlag+1;
        
    % Compute some statistics
    sige=sqrt((resid'*resid)/(nobs-nvarx));
    sigma2hat=sige(:).^2;
    invx=diag(inv(xt'*xt));
    bstd=sqrt(sigma2hat*invx);            

    betat=beta./bstd;
    pval=2*tcdf(-abs(betat),nobs-nvarx);                

    m0=eye(nobs)-1/nobs.*ones(nobs);   
    ee=diag(resid'*resid);                        
    R2=1-ee/(yt'*m0*yt);                        

    datest=dates(stIdx:enIdx);
    
    %% Then estimate the model over rolling sample
    T=size(datest,1);
    numRollingSamples=T-windowSize+1;

    cnt=0;
    betaB=zeros(numRollingSamples,nvarx);
    for t=1:numRollingSamples    
        [ytt,xtt] = constructArModel(y(stIdx:windowSize+stIdx+cnt-1),nlag,1,windowSize+cnt);
        [~,betat,~] = estimateArModel(ytt,xtt);    
        betaB(t,:)=betat;    
        cnt=cnt+1;
    end;

    %% Make figures of estimation results

    numa=20;
    residm=resid-repmat(mean(resid),[nobs 1]);
    [rho,~,ul,ll]=sacf(residm,numa,1);    
    f=figure('name','Autocorrelation in residuals');
    h=bar(rho,'k');
    set(gca,'xtick',1:numa);        
    hold on
    h1=plot(1:numa,ul,'color',[0.5 0.5 0.5],'linestyle','--');
    h2=plot(1:numa,ll,'color',[0.5 0.5 0.5],'linestyle','--');
    axis tight
    set(gca,'ylim',[-1 1]);

    ytm=yt-repmat(mean(yt),[nobs 1]);
    [rho,~,ul,ll]=sacf(ytm,numa,1);    
    f=figure('name','Autocorrelation in actual');
    h=bar(rho,'k');
    set(gca,'xtick',1:numa);        
    hold on
    h1=plot(1:numa,ul,'color',[0.5 0.5 0.5],'linestyle','--');
    h2=plot(1:numa,ll,'color',[0.5 0.5 0.5],'linestyle','--');
    axis tight
    set(gca,'ylim',[-1 1]);    

    f=figure('name','Fitted and actual');
    h=plot(yt,'k','linewidth',2);
    hold on
    h1=plot(yhat,'color',[0.5 0.5 0.5],'linewidth',2,'linestyle','--');                                                       
    if i==1
        lh=legend([h h1],{'GDP (YoY)','Fitted values'});
    else
        lh=legend([h h1],{'Inflation (YoY)','Fitted values'});
    end;
    set(lh,'box','off','location','northeast');    

    f=figure('name','Parameter stability');
    h=plot(betaB(:,1),'k','linewidth',2);
    hold on
    h1=plot(betaB(:,2),'color',[0.5 0.5 0.5],'linewidth',2);
    h2=plot(betaB(:,3),'k--','linewidth',2);
    lh=legend([h h1 h2],{'$\alpha$','$\phi_1$','$\phi_2$'});
    set(lh,'box','off','location','northwest','interpreter','latex');


end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yhat,betas,errors]=estimateArModel(yt,xt)    

% Simple OLS in Matlab
betas=xt\yt;
yhat=xt*betas;
errors=yt-yhat;                

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
function [yt,xt]=constructArModel(y,p,stIdx,enIdx)    

init=0;    
[nobs,nvar] = size(y);
xt=ones(nobs,nvar*p)*init;
for ii=1:p
    xt(1+ii:nobs,(nvar*(ii-1)+1):nvar*ii)=y(1:nobs-ii,:);
end    

yt=y(stIdx:enIdx);
xt=[ones(enIdx-stIdx+1,1) xt(stIdx:enIdx,:)];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bic=BIC(errors,p)

T=size(errors,1);
bic=log((errors'*errors)/T) + (p+1)*(log(T)/T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function aic=AIC(errors,p)

T=size(errors,1);
aic=log((errors'*errors)/T) + (p+1)*(2/T);
