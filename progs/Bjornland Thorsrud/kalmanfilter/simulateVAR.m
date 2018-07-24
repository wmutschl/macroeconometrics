function [Y,phi,sigma2,mu,nlag]=simulateVAR()
% PURPOSE: Simulate VAR data, estimate VAR, forecast and do IRF and VDC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s=RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);

%% Settings VAR(1)
nlag=1;
nvary=2;
T=101;
sigma2=[1 0.2;0.2 1];
%sigma2=[1 0;0 1];
mu=[1;1];
phi=[0.5 0;1 0.2];
%phi=[0.5 0;0 0.2];
y0=[2;1];
% random erros (draw from multivariate normal using chol decomp)
e=chol(sigma2,'lower')*randn(nvary,T);

y=[y0 nan(nvary,T-1)];
% generate VAR
for t=2:T    
    y(:,t)=mu+phi*y(:,t-1)+e(:,t);
end;
Y=y';