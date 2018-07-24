clear; 

% Quarterly data from FRED
% Gross National Product: Chain-type Price Index: Index 2005=100: SA, Q
% Sample 1954.Q4-2007.Q4
load gnpdeflator.txt; 
y = log(gnpdeflator(:,3));
T=size(y,1);
infl=(y(2:T,:)-y(1:T-1,:))*100;

pmax=12;
[sichat,hqchat,aichat]=pfind(infl,pmax);
p1 = aichat
p2 = 1;
numAutoCorr=p1+10;
OLSARp = ARpOLS(infl,p1,0,numAutoCorr);
OLSAR1 = ARpOLS(infl,p2,0,numAutoCorr);
[OLSARp.Qp OLSAR1.Qp]
r_k = nan(1,numAutoCorr);  % Initialize output vector
c0 = 1/(size(uhat,1))*(uhat' *uhat); % Compute variance
for h=1:numAutoCorr
    r_k(1,h) = 1/((size(uhat,1)-h)*c0) * (uhat(1+h:size(uhat,1),:)' * uhat(1:size(uhat,1)-h,:));
end
Q = size(uhat,1)*sum(r_k.*r_k);
Qp = chi2pdf(Q,h-p);