function [result,beta]=chapter9_consumptionIncomeCointegration()
%% PURPOSE: Engle-Granger two-step cointegration approach
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

%% Load data

data=xlsread('.\consumptionAndIncomeUS.xlsx','q');

pce=log(data(2,:))';
dspi=log(data(3,:))';
dates=data(1,:)';

%% Plot data
figure('name','Consumption and Income');
p1=plot(pce,'k','lineWidth',2);
hold on
p2=plot(dspi,'k--','lineWidth',2);                           
legend([p1 p2],{'Consumption','Income'},'location','northwest');
legend('boxoff');

%% Do cointegration regression
beta=[ones(size(dspi,1),1) dspi]\pce;
resid=pce-[ones(size(dspi,1),1) dspi]*beta;

figure('name','Consumption and Income residual');
p1=plot(resid,'k','lineWidth',2);

%% Do ADF test on data and estimated residual. 
data=[pce dspi resid];
vname={'pce','dspi','resid'};

p=-1:1; % no constant trend, constant, constant and trend
nlag=[2 4 6]; % ADF lags

nump=numel(p);
numnlag=numel(nlag);
for h=1:size(data,2)
    for i=1:nump
        for j=1:numnlag
            result(i,j,h)=unitRootTest(data(:,h),p(i),nlag(j),'vname',vname{h});            
        end;
    end;
end;

