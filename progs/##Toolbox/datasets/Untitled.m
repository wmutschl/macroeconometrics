Weberdata = readtable('WEBER.DAT','ReadVariableNames',0,'Format','%s   %f');
WeberDAT = table2cell(Weberdata)
xlswrite('WEBER.xlsx',WeberDAT)

[data,datatxt] = xlsread('WEBER.xlsx','Quarterly');
nobs = size(data,1);
% U:
% 1.) The monthly unemployment rate series from the Weber data is
% 	converted to quarterly frequency by averaging the monthly values.
% 2.) To account for a linear trend, the quarterly data are regressed
%     on a constant and a deterministic trend (using observations for 1948Q2-1987Q4).
% 3.) U is given by residuals from the regression in 2.).

y = data(:,2);
X = [ones(nobs,1) transpose(1:nobs)];
betta = inv(X'*X)*X'*y;
u = y-X*betta;

% DQ: 
% 1.) 1st differences of 100*log(real GNP) is taken, where "real GNP" is from the Weber data.
% 2.) The 1st differences are regressed on a constant and step dummy (1 after 1973Q4, 0 = else) 
%     in order to demean and adjust for the change in the output growth rate.
% 3.) DQ is given by residuals from the regression in 2.).
Dummy = zeros(nobs-1,1);
Dummy(find(strcmp(datatxt(2:end,1),'1973Q4'))+1:end)=1;
X = [ones(nobs-1,1) Dummy];
y = diff(100*log(data(:,1)));
betta = inv(X'*X)*X'*y;
u = y-X*betta;
