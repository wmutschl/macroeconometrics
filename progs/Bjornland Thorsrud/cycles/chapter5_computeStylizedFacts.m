function result=chapter5_computeStylizedFacts()
%% PURPOSE: Compute and compare different trend and cycles across different 
% time series
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
lambda=1600;
% Number of lags to be used to compute BN filter
nlag=8;
% Band pass filter: 
% pl - minimum period of oscillation of desired component 
% pu - maximum period of oscillation of desired component (2<=pl<pu<infinity).
pl=6;
pu=32;

ylim=[-1 1]; % for plotting cycles


%% Load data
ynames={'GDP','CONS','EXP','IMP','PROD','INVEST','EMP'};
[data,~]=xlsread('./stylizedFactsData.xlsx','q');

dates=data(1,:)';
y=data(2:end,:)';
yd=y(5:end,:)-y(1:end-4,:);

%% Do unit root testing

for h=1:size(y,2)   
    urTest(h)=unitRootTest(y(:,h),1,4,'vname',ynames{h});
end;

%% Compute Cycles and correlations
[T,N]=size(y);

[yc_hp,yc_li,yc_bp,yc_bn]=deal(nan(T,N));
for i=1:N
    [yc_hp(:,i),yc_li(:,i),yc_bp(:,i),yc_bn(:,i)]=computeCycle(y(:,i),lambda,nlag,pl,pu,T);
end;

%% Compute correlations on lead and lags

filterNameList={'li';'hp';'bp';'bn';'yd'};
ynamesLong=[];

leadAndLag=-4:1:4;

corrStylizedFacts=nan(N*N*4,numel(leadAndLag));
cnt=0;    
for i=1:N    
    for j=1:N
        [c_hp,~]=corrlela([yc_hp(:,i) yc_hp(:,j)],max(leadAndLag));
        [c_li,~]=corrlela([yc_li(:,i) yc_li(:,j)],max(leadAndLag));
        [c_bp,~]=corrlela([yc_bp(:,i) yc_bp(:,j)],max(leadAndLag));
        [c_bn,~]=corrlela([yc_bn(:,i) yc_bn(:,j)],max(leadAndLag));        
        [c_yd,~]=corrlela([yd(:,i) yd(:,j)],max(leadAndLag));        
        
        corrStylizedFacts(1+cnt:cnt+5,:)=[c_li';c_hp';c_bp';c_bn';c_yd'];             
        
        %filerNameList=cat(1,filerNameList,filerNameList);
        ynamesLong=cat(1,ynamesLong,strcat(ynames{i},'_',ynames{j},'_',filterNameList));
                
        cnt=cnt+5;                
    end;
end;

%% Make some correlation plots

f=figure('name','Stylized facts');

subplot(2,3,1)
idx=strcmpi('GDP_PROD_yd',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_PROD_yd','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

subplot(2,3,2)
idx=strcmpi('GDP_EMP_yd',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_EMP_yd','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

subplot(2,3,3)
idx=strcmpi('GDP_IMP_yd',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_IMP_yd','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

subplot(2,3,4)
idx=strcmpi('GDP_PROD_bn',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_PROD_bn','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

subplot(2,3,5)
idx=strcmpi('GDP_EMP_bn',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_EMP_bn','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

subplot(2,3,6)
idx=strcmpi('GDP_IMP_bn',ynamesLong);
h=bar(corrStylizedFacts(idx,:),'k');
title('GDP_IMP_bn','interpreter','none');
set(gca,'xticklabel',leadAndLag);
set(gca,'ylim',ylim);

%% result struct()

result.corrStylizedFacts=corrStylizedFacts;
result.urTest=urTest;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yc_hp,yc_li,yc_bp,yc_bn]=computeCycle(y,lambda,nlag,pl,pu,T)
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



