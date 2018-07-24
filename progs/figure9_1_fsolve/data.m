% DATA.M

% Quarterly data from FRED
% GDP: Chain-type Price Index: Index 2009
% Real Gross Domestic Product: Billions of Chained 2009 Dollars
% 1973.I-2013.II
load gdpdeflator2.txt; infl=dif(log(gdpdeflator2(56:end,3)))*100;
load realgdp2.txt; drgdp=dif(log(realgdp2(56:end,3)))*100;

% Monthly WTI spot price from Economagic, 1973.1-2007.12
load poil.txt; poilm=log(poil(:,3));
drpoil=(poilm(4:3:end)-poilm(1:3:end-3))*100-infl;
y=[drpoil infl drgdp];
