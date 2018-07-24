function chapter2_monteCarloArBias()
%% PURPOSE: Monte Carlo simulation of (OLS) AR. Estimation bias in small samples
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
phi = 0.9;
sigma2 = 0.5;

N = 500;
T = 10:1:250;
burn = 20;

%% Simulation
phiMCmean=mcsim(mu,phi,sigma2,N,T,burn);

%% Do the same, but now with:
phi2 = 0.01;
phiMCmean2=mcsim(mu,phi2,sigma2,N,T,burn);

%% Do the same, but now with:
phi3 = 1;
mu = 0;
phiMCmean3=mcsim(mu,phi3,sigma2,N,T,burn);

%% Make figures 
figure('name','MC simulation 1');
p1 = plot(phiMCmean,'k','lineWidth',2);
hold on
p2 = plot(phi.*ones(1,numel(T)+T(1)),'k--','lineWidth',2);
legend([p1 p2],{'Estimation','True'},'location','southeast');
legend('boxoff');
set(gca,'ylim',[0.4 1]);
set(gca,'ytick',0.4:0.1:1);
set(gca,'xlim',[0 260]);

figure('name','MC simulation 2');
p1 = plot(phiMCmean2,'k','lineWidth',2);
hold on
p2 = plot(phi2.*ones(1,numel(T)+T(1)),'k--','lineWidth',2);
legend([p1 p2],{'Estimation','True'},'location','southeast');
legend('boxoff');
set(gca,'ylim',[-0.4 0.2]);
set(gca,'ytick',-0.4:0.1:0.2);
set(gca,'xlim',[0 260]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phiMCmean=mcsim(mu,phi,sigma2,N,T,burn)

% Empty output
numT = numel(T);
phiMC = nan(N,numT+T(1)-1);

% Initial value - equal to unconditional mean of process
%y0 = mu/(1-phi);
y0 = 0;

for tT=T    
%     % the length of simulation sample for this T
%     tt = T(tT);    
    for s=1:N
        % simulate the process
        ys = nan(tT+burn,1);                    
        ys(1) = y0;        
        for t=2:(tT+burn)
            ys(t) = mu + phi*ys(t-1) + sqrt(sigma2)*randn(1,1);            
        end;        
        % Drop the burn-in observations 
        ys = ys(end-tT:end);        
        % Estimate coefficients
        ylag = ys(1:end-1);
        y = ys(2:end);
        const = ones(tT,1);
        % Easy OLS in Matlab
        coeffs = [const ylag]\y;                                
        % Save        
        phiMC(s,tT) = coeffs(2);        
    end;
end;
% Take the average over MC samples 
phiMCmean=mean(phiMC,1);








