% TABLE2_1.M

clear; 

% Quarterly data from FRED
% Gross National Product: Chain-type Price Index: Index 2005=100: SA, Q
% Real Gross Domestic Product: Billions of Chained 2005 Dollars: SAAR, Q
% Effective federal funds rate, M
% 1954.IV-2007.IV (exclude unconventional monetary policy after 2007.Q4)
load gnpdeflator.txt; infl=dif(log(gnpdeflator(:,3)))*100;
load realgnp.txt; drgdp=dif(log(realgnp(:,3)))*100;
load fedfunds.txt;
irate=[];
for i=1:3:length(fedfunds(:,3))
  irate=[irate; mean(fedfunds(i:i+2,3))];
end;    
xlswrite('test.xlsx',y)
% Same data set as in Rubio-Ramirez et al. (2010)
y=[drgdp irate infl]; [t,K]=size(y); 

% Find phat by AIC, HQC, and SIC for p=0,1,...,pmax
pmax=4;
[sichat,hqchat,aichat]=pfind(y,pmax)