function chapter5_computeTrendAndCycle()
%% PURPOSE: Compute and compare different trend and cycles
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

% HP filter smoothing paramter (3 alternatives)
lambda=[100 1600 40000];
% Number of lags to be used to compute BN filter
nlag=8;
% Band pass filter: 
% pl - minimum period of oscillation of desired component 
% pu - maximum period of oscillation of desired component (2<=pl<pu<infinity).
pl=6;
pu=32;

ylim=[-0.1 0.1]; % for plotting cycles

%% Load data
data=xlsread('./gdp.xlsx','q')';
y=log(data(:,2));
dates=data(:,1);
T=size(y,1);

%% Compute trend and cycles

% HP
[ytr_hp,yc_hp]=deal(nan(T,numel(lambda)));
for i=1:numel(lambda)
    ytr_hp(:,i)=hpfilter(y,lambda(i));
    yc_hp(:,i)=y-ytr_hp(:,i);
end;

% BN
[yc_bn,ytr_bn]=bnDecomp(y,nlag);

% BP
yc_bp=bpass(y,pl,pu);
ytr_bp=y-yc_bp;

% Linear
yc_li=detrend(y);
ytr_li=y-yc_li;

%% Plot output

f=figure('name','All: Cycle');
h=plot(zeros(T,1),'k','linewidth',2);
hold on
h1=plot(yc_li,'linewidth',2,'color',[0.5 0.5 0.5],'lineStyle','--');
h2=plot(yc_hp(:,2),'linewidth',2,'color',[0.5 0.5 0.5],'lineStyle','-');
h3=plot(yc_bn,'linewidth',2,'color',[0 0 0],'lineStyle','-');
h4=plot(yc_bp,'linewidth',2,'color',[0 0 0],'lineStyle','--');
set(gca,'ylim',ylim);
hl=legend([h1 h2 h3 h4],{'LT','HP-1600','BN','BP'});
set(hl,'box','off','location','northwest')

