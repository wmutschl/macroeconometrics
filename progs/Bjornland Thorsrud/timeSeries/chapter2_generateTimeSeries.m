function chapter2_generateTimeSeries()
%% PURPOSE: Generate different time series processes
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

%% Set the randstream (so that we get the same random numbers every time)
s=RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

%% Settings

% Number of autocorrelation lags to compute
numAutoCorr=21;


%% Generate AR
phi=[-0.8 0.4 0.9];
sigma=1;
mu=0;
T=200;
% Generate
y=nan(T,3);
y(1,1)=mu/(1-phi(1));
y(1,2)=mu/(1-phi(2));
y(1,3)=mu/(1-phi(3));
for t=2:T
    y(t,1)=mu+phi(1)*y(t-1,1)+randn()*sigma;    
    y(t,2)=mu+phi(2)*y(t-1,2)+randn()*sigma;        
    y(t,3)=mu+phi(3)*y(t-1,3)+randn()*sigma;            
end;

nameAr={'arn08phi','ar04phi','ar09phi'};
nameArAc={'arn08phiAc','ar04phiAc','ar09phiAc'};

for i=1:3
    plotseries(y(:,i),nameAr{i});
    plotautocorrTheory('ar',nameArAc{i},phi(i),numAutoCorr);
    plotautocorr(y(:,i),numAutoCorr,nameArAc{i});
end;

%% Generate MA
theta=[-0.5 0.5 0.9];
sigma=1;
mu=0;
T=200;
% Generate
y=nan(T,3);
y(1,1)=mu;
y(1,2)=mu;
y(1,3)=mu;

e=randn(T,3).*sigma;

for t=2:T
    y(t,1)=mu+e(t,1)+theta(1)*e(t-1,1);
    y(t,2)=mu+e(t,2)+theta(2)*e(t-1,2);
    y(t,3)=mu+e(t,3)+theta(3)*e(t-1,3);
end;

nameMa={'man05phi','ma05phi','ma09phi'};
nameMaAc={'man05phiAc','ma05phiAc','ma09phiAc'};

for i=1:3
    plotseries(y(:,i),nameMa{i});
    plotautocorrTheory('ma',nameMaAc{i},theta(i),numAutoCorr);
    plotautocorr(y(:,i),numAutoCorr,nameMaAc{i});
end;

%% Generate RW
sigma=1;
T=200;

% Generate
y=nan(T,3);y(1,:)=zeros(1,3);
for t=2:T
    y(t,1)=y(t-1,1)+randn()*sigma;    
    y(t,2)=y(t-1,2)+randn()*sigma;        
    y(t,3)=y(t-1,3)+randn()*sigma;            
end;
plotseries(y,'rw');
plotautocorr(y(:,1),numAutoCorr,'rwAc');

%% Generate RW with drift
sigma=1;
T=200;
mu=[-0.2 -0.01 0.2];

% Generate
y=nan(T,3);y(1,:)=zeros(1,3);
for t=2:T
    y(t,1)=mu(1)+y(t-1,1)+randn()*sigma;    
    y(t,2)=mu(2)+y(t-1,2)+randn()*sigma;        
    y(t,3)=mu(3)+y(t-1,3)+randn()*sigma;            
end;
plotseries(y,'rwDrift');
plotautocorr(y(:,1),numAutoCorr,'rwDriftAc');

%% Generate white noise
sigma=1;
T=200;
y=randn(T,1).*sigma;
plotseries(y,'whiteNoise');
plotautocorr(y,numAutoCorr,'whiteNoiseAc');

%% Generate Moving average
windowSize=5;
% Note, the filter command is Matlab code
y=filter(ones(1,windowSize)/windowSize,1,y);

plotseries(y,'whiteNoiseMovA');
plotautocorr(y,numAutoCorr,'whiteNoiseMovAAc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotseries(y,name)

f=figure('name','Time series plot');
if size(y,2)==3
    plot(y(:,1),'k','linewidth',2);
    hold on
    plot(y(:,2),'k--','linewidth',2);
    plot(y(:,3),'color',[0.5 0.5 0.5],'linewidth',2);
else
    h=plot(y,'k','linewidth',2);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotautocorrTheory(model,name,param1,hor)

acf=zeros(1,hor);
acf(1)=1;
switch model
    case {'ma'}
        acf(2)=param1/(1+param1^2);
    case {'ar'}
        for i=2:hor
            acf(i)=param1^(i-1);
        end;
end;

f=figure('name','Autocorrelation Theory');
h=bar(acf,'k');
%xlabel('k-values','fontsize',fontSize);    
%ylabel('sacf values','fontsize',fontSize);    
%title(ynames{i},'fontsize',fontSize);    
%legend('boxoff');    
set(gca,'xtick',1:hor);        
axis tight
set(gca,'ylim',[-1 1]);
set(gca,'xticklabel',num2str(str2num(get(gca,'xticklabel'))-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotautocorr(y,m,name)

nobs=size(y,1);
yy=y-repmat(mean(y),[nobs 1]);
[rho,~,ul,ll]=sacf(yy,m,1);    

f=figure('name','Autocorrelation');
h=bar(rho,'k');
set(gca,'xtick',1:m);        
hold on
h1=plot(1:m,ul,'color',[0.5 0.5 0.5],'linestyle','--');
h2=plot(1:m,ll,'color',[0.5 0.5 0.5],'linestyle','--');
axis tight
set(gca,'ylim',[-1 1]);
