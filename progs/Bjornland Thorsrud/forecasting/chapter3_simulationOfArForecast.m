function chapter3_simulationOfArForecast
%% PURPOSE: Simulate forecasts from AR processes. Generate density 
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
mu = 0.5;
phi = 0.8;
sigma2 = 0.10;

mu2 = 0.4;
phi1 =0.5;
phi2 = 0.2;

% The history
yt = [2.5 2.2 1.8 2.0 2.1 1.5 1.8 2]';

% Define the sign level, horizon and number of draws
percentiles = [0.10; 0.15; 0.25; 0.35]; % Same as Norges Bank report 
H = 10;
Nbig = 10000;
Nsmall = 100;

%% Simulate forecast 
ytF = forecastArModel([mu phi]',yt,1,H);
ytF2 = forecastArModel([mu2 phi1 phi2]',yt,2,H);
%% Generate uncertainty around the forecast

% First get the MSE...
mse = nan(H,1);
mse(1) = sigma2;
for h=2:H
    mse(h) = mse(h-1)+phi^(2*(h-1))*sigma2;    
end;
% ...then generate forecast intervals 
fanChart = nan(H,numel(percentiles)*2);
for pp=1:numel(percentiles);
    z2 = abs(norminv(percentiles(pp)/2,0,1));
    forcInt = nan(H,2);
    for h=1:H
        forcInt(h,1) = ytF(h)-z2*sqrt(mse(h));        
        forcInt(h,2) = ytF(h)+z2*sqrt(mse(h));
        
        fanChart(h,[pp (numel(percentiles)*2)-pp+1])=[forcInt(h,1) forcInt(h,2)];
    end;    
end;

nobs = numel(yt);

% Get quantile positions for plotting
quantLocationsDiffs=getQuantilePos(fanChart,yt,nobs);

%% Generate density forecast
forcDensityBig = nan(H,Nbig);
forcDensitySmall = nan(H,Nsmall);
for h=1:H
    forcDensityBig(h,:) = normrnd(ytF(h),sqrt(mse(h)),1,Nbig); 
    forcDensitySmall(h,:) = normrnd(ytF(h),sqrt(mse(h)),1,Nsmall); 
end;

%% Make output figures

% Make plot of point forecast
figure('name','Point forecast');
p1 = plot([yt;nan(H,1)],'k','lineWidth',2);
hold on
p2 = plot([nan(nobs-1,1);yt(end);ytF],'k--','lineWidth',2);
p3 = plot([nan(nobs-1,1);yt(end);ytF2],'k-*','lineWidth',2);
legend([p1 p2 p3],{'Actual','AR(1)','AR(2)'},'location','northwest');
legend('boxoff');
set(gca,'ylim',[1 4]);

% Make fanchart
f=figure('name','Fan chart');
fanchart=area(quantLocationsDiffs);
%set colormap
setQuantileCol(fanchart);
hold on 
p1 = plot([yt;nan(H,1)],'k','lineWidth',2);
p2 = plot([nan(nobs-1,1);yt(end);ytF],'k--','lineWidth',2);
legend([p1 p2],{'Actual','Point forecast and uncertainty'},'location','northwest');
legend('boxoff');
set(gca,'ylim',[1 4]);

% Make histograms of density along with normal curve
figure('name','Histogram with along with normal: Snmall N')
h = histfit(forcDensitySmall(1,:)');
set(h(1),'FaceColor',[0.7 0.7 0.7]);
set(h(2),'Color',[0 0 0]);

figure('name','Histogram with along with normal: Large N')
h = histfit(forcDensityBig(1,:)');
set(h(1),'FaceColor',[0.7 0.7 0.7]);
set(h(2),'Color',[0 0 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yhatF=forecastArModel(betas,yt,p,H)

% initial conditioning vector 
ycond = yt(end:-1:end-p+1)';

yhatF = nan(H,1);
for h=1:H    
    % add the constant to the conditioning information and forecast
    yhatF(h) = [1 ycond]*betas;    
    % update the conditioning vector with previous periods forecast
    ycond = [yhatF(h) ycond(1:end-1)];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

