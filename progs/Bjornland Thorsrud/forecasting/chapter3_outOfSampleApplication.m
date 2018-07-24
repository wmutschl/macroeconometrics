function chapter3_outOfSampleApplication
%% PURPOSE: Run an out-of-sample forecasting experiment using Norwegian 
% GDP growth. 
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
%
% 1) Define the settings of the experiment
% 2) Load data and transform
% 3) DO AIC and BIC test for lag length
% 4) Do out-of-sample experiment
% 5) Evalaute individual models and combination
% 6) Make figures and tables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1) Settings

% time index settings
oosStart = 1990.01; % start of out-of-sample forecasting
oosEnd = 2013.01;   % end of out-of-sample forecasting
eStart = 1980.01;   % first period for estimation 
oosEvalStart = 1993.01; % start date for evaluation of the out-of-sample exp.
oosEvalEnd = 2013.01;   % end date for evaluation of the out-of-sample exp.

% Maximum number of lags to consider in model
P = 4;
% Maximum nuber of forecasts (horizon)
H = 4;

% Choose to do rolling or not rolling estimation of AR models in combo.
% Equals false if expanding window should be used
doRollingEstimation = true;

% Set the number of lags to be used for the two AR models
arModelLagInCombo = [1 4]; 

%% 2) Load data. Construct growth rate

[data,~] = xlsread('./gdpNorway.xlsx','gdp');
dates = data(:,1);
yLevel = data(:,2);
% take log and get the log year on year diff (scaled by 100 to interpret as
% percent)
yLevelLog = log(yLevel).*100;
y = yLevelLog(5:end,1)-yLevelLog(1:end-4,1);
%y=diff(yLevelLog);
% correct dates for one lost observation due to diff
dates = dates(5:end);

% Do some simple error checking on the settings 
if isempty(find(dates==eStart))
    error('Start of estimation sample not in loaded dates')
end;
if find(dates==eStart)-P<1
    error('P is too large relative to eStart')    
end;
if find(dates==eStart)-max(arModelLagInCombo)<1
    error('Max of arModelLagInCombo is too large relative to eStart')    
end;

%% 3) Do aic and bic test to choose lag length. Do this on whole sample and on sample
%% up to oosStart; which is the information we would actually be having in oosStart
[aicTest,bicTest] = deal(nan(P,2)); % empty output variables

% For oosStart
stIdx = find(dates==eStart);        % start index for estiamtion
enIdx = find(dates==oosStart);      % end index for estimation
for p=1:P    
    [yt,xt] = constructArModel(y,p,stIdx,enIdx);
    [~,~,errors] = estimateArModel(yt,xt);    
    bicTest(p,1) = BIC(errors,p);
    aicTest(p,1) = AIC(errors,p);        
end;
% For whole sample
stIdx = find(dates==eStart);        % start index for estiamtion
enIdx = find(dates==oosEnd);      % end index for estimation
for p=1:P    
    [yt,xt] = constructArModel(y,p,stIdx,enIdx);
    [~,~,errors] = estimateArModel(yt,xt);    
    bicTest(p,2) = BIC(errors,p);
    aicTest(p,2) = AIC(errors,p);        
end;

%% 4) Do out-of-sample forecasting experiment. We use one AR(1) model, and 
%% one AR(4) model as well as a combination of the two 

% add H nan values to y such that the evaluation do not run out of
% "observations"
yLong = [y;nan(H,1)];

% make some time indexes 
oosStartIdx = find(dates==oosStart);      
oosEndIdx = find(dates==oosEnd);      
eStartIdx = find(dates==eStart);       

tildeT = oosEndIdx-oosStartIdx+1;
T = oosStartIdx;
barT = size(dates,1);

% empty output
[errors,forecasts] = deal(nan(barT,3,H)); % for individual models and combo
[mse,weights] = deal(nan(barT,2,H));     % for individual models only

% Do the experiment for the individual models
for starT=T:barT    
    % Realizations for forecasts given at time starT
    yObs = yLong(starT+1:starT+H);
    
    % AR model 1
    if ~doRollingEstimation
        % Expanding estimation window
        [yt,xt] = constructArModel(y,arModelLagInCombo(1),eStartIdx,starT);    
    else
        % Roling estimation window
        [yt,xt] = constructArModel(y,arModelLagInCombo(1),eStartIdx+(starT-T),starT);    
    end;
    
    [~,betas,~] = estimateArModel(yt,xt);            
    yhatF1 = forecastArModel(betas,yt,arModelLagInCombo(1),H);              
    
    % output
    forecasts(starT,1,:) = yhatF1;
    errors(starT,1,:) = yObs-yhatF1;                
    
    % AR model 2
    if ~doRollingEstimation
        % Expanding estimation window    
        [yt,xt] = constructArModel(y,arModelLagInCombo(2),eStartIdx,starT);
    else
        % Roling estimation window
        [yt,xt] = constructArModel(y,arModelLagInCombo(2),eStartIdx+(starT-T),starT);        
    end;
    
    [~,betas,~] = estimateArModel(yt,xt);    
    yhatF4 = forecastArModel(betas,yt,arModelLagInCombo(2),H);                
    
    % output
    forecasts(starT,2,:) = yhatF4;
    errors(starT,2,:) = yObs-yhatF4;                                                             
end;

% Compute mse based on the errors, and get weights
for h=1:H
    mse(T:barT-h,1,h) = (cumsum(errors(T:barT-h,1,h).^2)./(1:tildeT-h)');
    mse(T:barT-h,2,h) = (cumsum(errors(T:barT-h,2,h).^2)./(1:tildeT-h)');    
    
    weights(T:barT-h,1,h) = mse(T:barT-h,1,h)./(mse(T:barT-h,1,h)+mse(T:barT-h,2,h));
    weights(T:barT-h,2,h) = mse(T:barT-h,2,h)./(mse(T:barT-h,1,h)+mse(T:barT-h,2,h));    
end;

% Do model combination
for h=1:H
    for starT=(T+h):barT    
        % Realizations for forecasts given at time starT
        yObs = yLong(starT+h);    
    
        % Linear combination
        yhatFCombo = weights(starT-h,1,h)*forecasts(starT,1,h) + weights(starT-h,2,h)*forecasts(starT,2,h);
        
        % Output
        forecasts(starT,3,h) = yhatFCombo;
        errors(starT,3,h) = yObs-yhatFCombo;            
    end;
end;

%% 5) Evaluate individual models and combo over same sample according to RMSE

% indexing
oosEvalStartIdx = find(dates==oosEvalStart);      
oosEvalEndIdx = find(dates==oosEvalEnd);      

% empty output
[rmse,bias] = deal(nan(H,3));
 
% Compute rmse for indiviudal models and combo
for h=1:H    
    bias(h,1) = mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,1,h));
    bias(h,2) = mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,2,h));
    bias(h,3) = mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,3,h)); 
    
    rmse(h,1) = (mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,1,h).^2))^0.5;
    rmse(h,2) = (mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,2,h).^2))^0.5;
    rmse(h,3) = (mean(errors(oosEvalStartIdx:oosEvalEndIdx-h,3,h).^2))^0.5;                
end;

%% 6) Make output figures and table 

% Make plot combined forecast along with actual 
figure('name','Forecast and actual');
p1=plot(y(oosEvalStartIdx+1:oosEvalEndIdx),'k','lineWidth',2);
hold on
p2=plot(forecasts(oosEvalStartIdx:oosEvalEndIdx-1,3,1),'k--','lineWidth',2);
legend([p1 p2],{'Actual','Forecast'});
legend('boxoff');

% Make plot of the weights
figure('name','Weights');
p1=plot(weights(oosStartIdx:oosEndIdx,1,1),'k','lineWidth',2);
hold on
p2=plot(weights(oosStartIdx:oosEndIdx,2,1),'k--','lineWidth',2);
legend([p1 p2],{'Model 1','Model 2'});
legend('boxoff');

% Make hairy plot of Model 1 
hndl=figure('name','Hairy plot');
for t=1:tildeT
    hairy=[y(oosStartIdx:oosStartIdx+t-1);permute(forecasts(oosStartIdx+t-1,1,:),[3 2 1]);nan(tildeT-t-H,1)];
    p2=plot(hairy,'k--','lineWidth',2);
    hold on
end;
p1=plot(y(oosStartIdx:oosEndIdx),'k','lineWidth',2);
legend([p1 p2],{'Actual','Model 1'});
legend('boxoff');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra function used by the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yhatF=forecastArModel(betas,yt,p,H)

% initial conditioning vector 
ycond=yt(end:-1:end-p+1)';

yhatF=nan(H,1);
for h=1:H    
    % add the constant to the conditioning information and forecast
    yhatF(h)=[1 ycond]*betas;    
    % update the conditioning vector with previous periods forecast
    ycond=[yhatF(h) ycond(1:end-1)];    
end
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

