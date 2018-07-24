function [result_op,result_cl]=chapter8_estimateVARmp()
% PURPOSE: Estimate standard monetary policy VAR for Sweden
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
nlag=3;
% Impulse-response horizon
nirf=24;

%% Run 
[result_op,result_cl]=computeModel(nlag,nirf);

%% Make irf figures for closed economy using normal ordering

irf=result_cl.irf;
% Scale MP shock to 1pp
normConstant=irf(3,3,1);
irf(:,3,:)=irf(:,3,:)./normConstant;
% make level of GDP
irf(1,1,:)=cumsum(irf(1,1,:),3);

makeIRFfigure(permute(irf(1,3,:).*100,[3 1 2]),[],{'GDP','MP'},'ch'); % .*100 to get \%
makeIRFfigure(permute(irf(2,3,:).*100,[3 1 2]),[],{'INF','MP'},'ch'); % .*100 to get \%
makeIRFfigure(permute(irf(3,3,:),[3 1 2]),[],{'MP','MP'},'ch');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_op,result_cl]=computeModel(nlag,nirf)


%% Load data
[data,txt]=xlsread('bjornland_swe_Data.xlsx','Sheet1','A5:P98');
txt=[txt(1:2,2:end) {'c';'c'} {'tr';'tr'}];                

T=size(data,1);
tr=(1:T)';
c=ones(T,1);
data=[data c tr];

% endogen variables
endogc={'GDP','Aldp','r'}; % closed economy
endogo={'rtwi','GDP','Aldp','r','RERinv'}; % open economy
% exogen variables
exog={'c','tr','du93Q1','du95Q4','du92Q3'};       
        
%% extract data

% for closed economy model
yc=nan(T,numel(endogc));
for i=1:numel(endogc)
    idx=strcmpi(endogc{i},txt(2,:));
    if i==1 
        yc(:,i)=log(data(:,idx));
    else
        yc(:,i)=data(:,idx);        
    end;
end;
% for open economy model
yo=nan(T,numel(endogo));
for i=1:numel(endogo)
    idx=strcmpi(endogo{i},txt(2,:));
    if i==2
        yo(:,i)=log(data(:,idx));
    elseif i==5
        yo(2:end,i)=diff(log(data(:,idx)));
        %yo(:,i)=log(data(:,idx));
    else
        yo(:,i)=data(:,idx);        
    end;
end;
e=nan(T,numel(exog));
for i=1:numel(exog)
    idx=strcmpi(exog{i},txt(2,:));
    e(:,i)=data(:,idx);
end;

%% Closed economy estimation. Cholesky identification
x_tmp=[latMlag(yc,nlag) e];
y=yc(2+nlag:end,:);
x=x_tmp(2+nlag:end,:);
nvary=size(y,2);

result_cl=estimateVAR(y,x,nvary,nlag);
% Point estimate mse
betac=varGetCompForm(result_cl.beta,[],result_cl.nlag,result_cl.nvary);
result_cl.mse=getMse(betac,result_cl.Sigma2,nirf);
% Cholesky identification
result_cl.irf=varImpulseResponsFast(betac,result_cl.nvary,result_cl.nlag,nirf,chol(result_cl.Sigma2,'lower'));
result_cl.vdc=getVdcFromIrf(result_cl.irf);
    
%% Open economy estimation. Cholesky and long run identiication
x_tmp=[latMlag(yo,nlag) e];
y=yo(2+nlag:end,:);
x=x_tmp(2+nlag:end,:);
nvary=size(y,2);

result_op=estimateVAR(y,x,nvary,nlag);
% Point estimate mse
betac=varGetCompForm(result_op.beta,[],result_op.nlag,result_op.nvary);
result_op.mse=getMse(betac,result_op.Sigma2,nirf);
% Cholesky identification
result_op.irf=varImpulseResponsFast(betac,result_op.nvary,result_op.nlag,nirf,chol(result_op.Sigma2,'lower'));
result_op.vdc=getVdcFromIrf(result_op.irf);










