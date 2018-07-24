function chapter1_monteCarloOLS()
%% PURPOSE: Monte Carlo simulation of OLS estimator
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
b0 = 0.5;
b1 = 1;
b2 = 0.2;
Sigma = [0.7 0.4 0;
       0.4 0.3 0;
       0  0  1];
T=[100 500 1000];
N=5000;

% Generate x values
sv = chol(Sigma(1:2,1:2))'*randn(2,max(T));    
x1 = sv(1,:)';
x2 = sv(2,:)';
x0 = ones(max(T),1);

%% Monte Carlo simulation of baseline model

% Empty output
paramSim=nan(N,3,numel(T));
paramSim2=nan(N,2,numel(T));

for t=1:numel(T)
    for s=1:N
        % generate random errors
        u = randn(T(t),1);    
        % simulate the process
        y = b0*x0(1:T(t)) + b1*x1(1:T(t)) + b2*x2(1:T(t)) + u;           
        % Estimate coefficients by OLS (easy in Matlab) and save
        paramSim(s,:,t) = [x0(1:T(t)) x1(1:T(t)) x2(1:T(t))]\y;                                
        % Mis-specified model         
        paramSim2(s,:,t) = [x0(1:T(t)) x1(1:T(t))]\y;                                
    end;
end;

%% Monte Carlo simulation alternative

% Empty output
paramSim3=nan(N,2,numel(T));
paramSim4=nan(N,3,numel(T));

for t=1:numel(T)
    for s=1:N
        % generate random errors
        u = randn(T(t),1);    
        % simulate the process
        y = b0*x0(1:T(t)) + b1*x1(1:T(t)) + u;           
        % Correct estimation
        paramSim3(s,:,t) = [x0(1:T(t)) x1(1:T(t))]\y;                                                
        % Mis-specified model         
        paramSim4(s,:,t) = [x0(1:T(t)) x1(1:T(t)) x2(1:T(t))]\y;                                        
    end;
end;


%% Make output figures

figure('name','Histogram of b1')
h = histfit(paramSim(:,2,1)');
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
hold on
h2 = histfit(paramSim(:,2,3)');
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');  
line([b1 b1],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

figure('name','Histogram of b2')
h = histfit(paramSim(:,3,1)');
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
hold on
h2 = histfit(paramSim(:,3,3)');
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');
line([b2 b2],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

figure('name','Histogram of b1 - bias')
h = histfit(paramSim2(:,2,1)');
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
hold on
h2 = histfit(paramSim2(:,2,3)');
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');
line([b1 b1],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

figure('name','Histogram of b1 - variance')
h = histfit(paramSim4(:,2,1)');
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'Color',[0 0 0],'lineStyle','--');
hold on
h2 = histfit(paramSim3(:,2,1)');
set(h2(1),'FaceColor',[0.95 0.95 0.95]);
set(h2(2),'Color',[0 0 0],'lineStyle','-');
line([b1 b1],[0 max(get(gca,'ylim'))],'color',[0 0 0],'linewidth',2)

