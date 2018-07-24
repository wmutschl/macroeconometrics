function chapter6_simultaneityMonteCarlo()
% PURPOSE: Do Monte Carlo simulation to show OLS versus TSLS when 
% simultaneity is a problem
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
nobs = 200;

x1 = randn(nobs,1);
x2 = randn(nobs,1);
b1 = 1.0;
b2 = 0.0;   %1.0
iota = ones(nobs,1);

y1 = zeros(nobs,1);
y2 = zeros(nobs,1);
evec = randn(nobs,1);

% create simultaneously determined variables y1,y2
for i=1:nobs;
    y1(i,1) = iota(i,1)*1.0 + x1(i,1)*b1 + evec(i,1);
    y2(i,1) = iota(i,1)*1.0 + y1(i,1)*1.0 + x2(i,1)*b2 + evec(i,1);
end;
          
%% Simple estimation

% use all exogenous in the system as instruments
xall = [iota x1 x2];            

% do ols regression
x=[y1 iota x2];
beta_ols=x\y2;

% do tsls regression
result=tslsEstimation(y2,y1,[iota x2],xall); 

disp('OLS and TSLS estimates:')
disp('OLS')
disp(beta_ols')
disp('TSLS')
disp(result.beta')

%% Monte Carlo

% do Monte Carlo looping
niter=1000;
N=[100 500 1000];

x1B=randn(max(N),1);
x2B=randn(max(N),1);
iotaB=ones(max(N),1);    
xallB=[iotaB x1B x2B];            

bolsB=cell(1,numel(N));
b2slsB=cell(1,numel(N));
tolsB=cell(1,numel(N));
t2slsB=cell(1,numel(N));

for d=1:numel(N)
       
    nobs=N(d);    
        
    bols=nan(size(beta_ols,1),niter);
    b2sls=nan(size(beta_ols,1),niter);    
    tols=nan(size(beta_ols,1),niter);
    t2sls=nan(size(beta_ols,1),niter);        
    
    for iter=1:niter;

        y1 = zeros(nobs,1);
        y2 = zeros(nobs,1);
        evec = randn(nobs,1);

        % create simultaneously determined variables y1,y2
        for i=1:nobs;
            y1(i,1) = iotaB(i,1)*1.0 + x1B(i,1)*b1 + evec(i,1);
            y2(i,1) = iotaB(i,1)*1.0 + y1(i,1)*1.0 + x2B(i,1)*b2 + evec(i,1);
        end;

        % do ols regression
        results_ols_mc=olsEstimation(y2,[y1 iotaB(1:nobs,:) x2B(1:nobs,:)]);

        % do tsls regression        
        results_tsls_mc=tslsEstimation(y2,y1,[iotaB(1:nobs,:) x2B(1:nobs,:)],xallB(1:nobs,:)); 

        bols(:,iter)=results_ols_mc.beta;
        b2sls(:,iter)=results_tsls_mc.beta;
        tols(:,iter)=results_ols_mc.tstat;
        t2sls(:,iter)=results_tsls_mc.tstat;        
        
    end;
    
    bolsB{d}=bols;
    b2slsB{d}=b2sls;
    tolsB{d}=tols;
    t2slsB{d}=t2sls;    
    
end;

%% Plot some convergence results

figure('name','Histogram of y1 coefficient: OLS')
h = histfit(bolsB{1}(1,:));
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
set(gca,'xlim',[0.9 1.7]);
hold on
h2 = histfit(bolsB{3}(1,:));
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');
line([1 1],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

figure('name','Histogram of y1 coefficient: TSLS')
h = histfit(b2slsB{1}(1,:));
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
hold on
h2 = histfit(b2slsB{3}(1,:));
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');
line([1 1],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

