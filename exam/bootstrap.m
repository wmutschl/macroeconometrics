clear; close all; clc
% Bootstrap der Residuen
% Generate data
TT=100;
u = exprnd(1,TT,1)-1;
y = nan(TT,1);
c = 1;
phi = 0.8;
y(1) = 1/(1-0.8);
for t=2:TT
    y(t) = c + phi*y(t-1) + u(t);
end
% OLS estimation
OLS = ARpOLS(y,1,1);
uhat = OLS.uhat;
chat = OLS.thetahat(1);
phihat = OLS.thetahat(2);
sig_phihat = OLS.sig_thetahat(2);
tstat = phihat/sig_phihat;

% Perzentil-t-Bootstrap
R = 10000;
taustar = nan(R,1);
for r=1:R
    %randi(10,1,10)
    %ustar = datasample(uhat,T_eff);
    ustar = uhat(ceil(size(uhat,1)*rand(OLS.T_eff,1)),:);
    ystar = nan(size(y,1),1);
    ystar(1) = y(1);     % Intialize the first p observations with real data
    for t=2:size(y,1)
        ystar(t) = chat+phihat*ystar(t-1)+ustar(t-1,1);
    end
    OLSstar = ARpOLS(ystar,1,1);
    phistar = OLSstar.thetahat(2);
    sig_phistar = OLSstar.sig_thetahat(2);
    taustar(r) = (phistar-phihat)/sig_phistar;
end

taustar = sort(taustar);
Lower_Boot = phihat-taustar(0.975*R)*sig_phihat;
Upper_Boot = phihat-taustar(0.025*R)*sig_phihat;
Lower_Approx = phihat-1.96*sig_phihat;
Upper_Approx = phihat+1.96*sig_phihat;
table([Lower_Approx;Lower_Boot],[Upper_Approx;Upper_Boot],'RowNames',{'Approx' 'Boot'},'VariableNames',{'Lower' 'Upper'})

