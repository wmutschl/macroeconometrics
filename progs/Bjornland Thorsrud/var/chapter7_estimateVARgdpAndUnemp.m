function result=chapter7_estimateVARgdpAndUnemp()
% PURPOSE: Estimate VAR using Blanchard and Quah data
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

% For lag selection test
maxLag=8;

% Forecast horizon and percentiles
H=10;
percentiles = [0.10; 0.15; 0.25; 0.35]; % Corresponding to Norges Bank report


%% Load data

[data,txt]=xlsread('./blanchardquaData.xlsx','data');

y1=data(:,2);
y2=data(:,3);
dates=data(:,1);

Y=[y1 y2];
[T,nvary]=size(Y);
ynames={'GDP','U'};


%% Estimate VAR by OLS and do lag selection test

% Do lag selection test
bicAndAic=nan(maxLag,2);
for i=1:maxLag
    Xt=[ones(T,1) latMlag(Y,i)];
    % use same number of observations for all estimations
    yt=Y(1+maxLag:end,:);
    xt=Xt(1+maxLag:end,:);

    nobst=size(yt,1);    
    % OLS
    coefft=(xt\yt)';
    residt=yt-(coefft*xt')';

    epe=(residt'*residt)/nobst;
    
    aicVal=aicTest(epe,nvary,i,nobst);
    bicVal=bicTest(epe,nvary,i,nobst);
    
    bicAndAic(i,:)=[bicVal aicVal];
end;
% Get number of lags based on lowest criteria...
[mi,li]=min(bicAndAic);
[mi,col]=min(mi);
nlag=find(bicAndAic(:,col)==mi);

%% Estimate by OLS and compute some statistics using lag length from above
[y,x,alfa,beta,yhat,resid,sigma2hat,bstd,R2,SS,nobs]=varEst(T,Y,nlag);
% Compute Granger causality
[gr_ftest,gr_fprob]=pgranger(y,x,true,nlag,nvary,resid);
% Put results into result struct
result.alfa=alfa;
result.beta=beta;
result.yhat=yhat;
result.y=y;
result.x=x;
result.resid=resid;
result.sigma2hat=sigma2hat;
result.bstd=bstd;
result.R2=R2;
result.gr_ftest=gr_ftest;
result.gr_fprob=gr_fprob;

%% Estimate by OLS and compute some statistics using 8 lags
[yl,xl,~,~,~,residl,~,~,R2l,~,~]=varEst(T,Y,8);
% Compute Granger causality
[gr_ftestl,gr_fprobl]=pgranger(yl,xl,true,8,nvary,residl);
% Put results into result struct
result.gr_ftestl=gr_ftestl;
result.gr_fprobl=gr_fprobl;

%% Generate forecasts
yhatF=forecastVARModel([alfa beta]',y,nlag,H,nvary);
% Put results into result struct
result.yhatF=yhatF;

% Compute MSE
[betac,~]=varGetCompForm(beta,alfa,nlag,nvary);
Sigma2=(resid'*resid)/nobs;
mse=getMse(betac,Sigma2,10);
% Put results into result struct
result.mse=mse;

% Generate uncertainty around the forecast
fanChart = nan(H,numel(percentiles)*2,nvary);
for j=1:nvary
    for pp=1:numel(percentiles);
        z2 = abs(norminv(percentiles(pp)/2,0,1));
        forcInt = nan(H,2);
        for h=1:H
            forcInt(h,1) = yhatF(h,j)-z2*sqrt(mse(j,j,h));        
            forcInt(h,2) = yhatF(h,j)+z2*sqrt(mse(j,j,h));

            fanChart(h,[pp (numel(percentiles)*2)-pp+1],j)=[forcInt(h,1) forcInt(h,2)];
        end;    
    end;
end;

%% Check the residuals for autocorrelation (informally by plotting)
residm=resid-repmat(mean(resid),[nobs 1]);
numa=20;
for i=1:nvary
    [rho,~,ul,ll]=sacf(residm(:,i),numa,1);    
    
    f=figure('name','Autocorrelation in residuals');
    h=bar(rho,'k');
    set(gca,'xtick',1:numa);        
    hold on
    h1=plot(1:numa,ul,'color',[0.5 0.5 0.5],'linestyle','--');
    h2=plot(1:numa,ll,'color',[0.5 0.5 0.5],'linestyle','--');
    axis tight
    set(gca,'ylim',[-1 1]);   
end;

%% Plot data

f=figure('name','Time series plot: GDP data');
h=plot(Y(1+nlag:end,1),'k','linewidth',2);
hold on
h1=plot(yhat(:,1),'color',[0.5 0.5 0.5],'linewidth',2,'linestyle','--');
set(gca,'ylim',[-4 4.5]);
hl=legend([h h1],{'GDP','yhat'});
set(hl,'box','off','location','northwest')

f=figure('name','Time series plot: U data');
h=plot(Y(1+nlag:end,2),'k','linewidth',2);
hold on
h1=plot(yhat(:,2),'color',[0.5 0.5 0.5],'linewidth',2,'linestyle','--');
set(gca,'ylim',[-4 4.5]);
hl=legend([h h1],{'Unemployment','yhat'});
set(hl,'box','off','location','northwest')

% Make fanchart
% Get quantile positions for plotting
nobs=8;
for j=1:nvary    
    yt=y(end-nobs+1:end,j);
    quantLocationsDiffs=getQuantilePos(fanChart(:,:,j),yt,nobs);

    figure('name','Fan chart');
    fanchart=area(quantLocationsDiffs);
    %set colormap
    setQuantileCol(fanchart);
    hold on 
    p1 = plot([yt;nan(H,1)],'k','lineWidth',2);
    p2 = plot([nan(nobs-1,1);yt(end);yhatF(:,j)],'k--','lineWidth',2);   
    legend([p1 p2],{'Actual','Point forecast and uncertainty'},'location','northwest');
    legend('boxoff');   
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,x,alfa,beta,yhat,resid,sigma2hat,bstd,R2,SS,nobs]=varEst(T,Y,nlag)

X=[ones(T,1) latMlag(Y,nlag)];
y=Y(1+nlag:end,:);
x=X(1+nlag:end,:);

[nobs,nvary]=size(y);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yhatF=forecastVARModel(betas,yt,p,H,nvary)

% initial conditioning vector 
ycond = yt(end:-1:end-p+1,:)';
ycond = ycond(:)';
yhatF = nan(H,nvary);

for h=1:H    
    % add the constant to the conditioning information and forecast
    yhatF(h,:) = [ones(1,1) ycond]*betas;    
    % update the conditioning vector with previous periods forecast
    ycond = [yhatF(h,:) ycond(:,1:end-nvary)];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quantLocationsDiffs=getQuantilePos(fanChart,yt,nobs)

quantLocations=fanChart';
quantLocationsDiffs=(quantLocations(2:end,:)-quantLocations(1:end-1,:)); % (h x q)
quantLocationsDiffs=[quantLocations(1,:);quantLocationsDiffs]';
% add obseravble
qc=size(quantLocationsDiffs,2);
quantLocationsDiffs(isnan(quantLocationsDiffs))=0;
quantLocationsDiffs=[[yt zeros(nobs,qc-1)];quantLocationsDiffs];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setQuantileCol(fanchart)

%colormap cool;
numAreas=length(fanchart);
st=ceil((numAreas-1)/2);
cmap=vivid(st,'k',[0.70 0.96]);
cmap=[flipud(cmap(1:st,:));cmap(2:st,:)];
set(fanchart(1),'FaceColor','none','linestyle','none','EdgeColor','none');	 
for i=2:numAreas
    set(fanchart(i),'FaceColor',cmap(i-1,:),'linestyle','none','EdgeColor','none')	 
end