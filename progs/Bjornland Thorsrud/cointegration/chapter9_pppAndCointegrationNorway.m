function [result,adf_real_s]=chapter9_pppAndCointegrationNorway()
%% PURPOSE: Run PPP model and check for cointegration in a VECM
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

%% Get data

[data,txt]=xlsread('./pppDataNorway.xlsx','data');

dates=data(:,1);
% All data are in logs 
lp=data(:,3);
lpf=data(:,4);
ls=data(:,2);
real_s=data(:,5);
% seasonal dummies
s1=data(:,6);
s2=data(:,7);
s3=data(:,8);

nvary=3;
vnames={'ls','lp','lpf'};

%% Make plots of data

f=figure('name','all');
h1=plot(lpf./lp,'k','linewidth',2);
hold on
h2=plot(ls,'color',[0.5 0.5 0.5],'linewidth',2);
h3=plot(real_s,'color',[0.5 0.5 0.5],'linewidth',2,'lineStyle','--');
lh=legend([h1 h2 h3],{'$CPI/CPI^{f}$','Nominal exchange rate','Real exchange rate'});
set(lh,'box','off','location','northeast','interpreter','latex');

%% Do ADF test

idx=find(dates==1997.04);

adf_lp=unitRootTest(lp(1:idx),1,7,'vname','lp');
adf_lpf=unitRootTest(lpf(1:idx),1,7,'vname','lpf');
adf_ls=unitRootTest(ls(1:idx),0,7,'vname','ls');
adf_real_s=unitRootTest(real_s(1:idx),0,7,'vname','real_s');

%% Estimate VAR in levels and check stationarity
nlag=5;

y=[ls lp lpf];
x=latMlag(y,nlag);
y=y(nlag+1:end,:);
T=size(y,1);
c=ones(T,1);
tr=(1:T)';
sd=[s1 s2 s3];
sd=sd(nlag+1:end,:);
x=[x(nlag+1:end,:) c tr sd];

% adjust for sample
datess=dates(nlag+1:end);
idxs=find(datess==1972.02);
idxe=find(datess==1997.04);
yi=y(idxs:idxe,:);
xi=x(idxs:idxe,:);

result=estimateVAR(yi,xi,nvary,nlag);

% get eigenvalues of comp form beta
betac=varGetCompForm(result.beta,[],nlag,nvary);

disp('Eigenvalues of comp form matrix')
disp(eig(betac))

%% Do my VECM (restricted)
y=[ls lp lpf];
x=latMlag(y,nlag);

yi=y(nlag+1:end,:);

yd=diff(yi);
xd=diff(x);

T=size(yd,1);
c=ones(T,1);
tr=(1:T)';

sd=[s1 s2 s3];
sd=sd(nlag+2:end,:);

vecm=y*[1 -1 1]'; % enforcing PPP to hold
vecm=latMlag(vecm,1);
vecm=vecm(nlag+2:end,:);

xd=[xd(nlag+1:end,:) c tr sd vecm];

% adjust for sample
datess=dates(nlag+2:end);
idxe=find(datess==1997.04);
idxs=find(datess==1972.02);
yi=yd(idxs:idxe,:);
xi=xd(idxs:idxe,:);

% output
result=estimateVAR(yi,xi,nvary,nlag);





