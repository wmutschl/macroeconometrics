% Copyright Andrew Binning 2013
% Please feel free to use and modify this code as you see if fit. If you
% use this code in any academic work, please cite 
% Andrew Binning, 2013.
% "Underidentified SVAR models: A framework for combining short and long-run restrictions with sign-restrictions,"
% Working Paper 2013/14, Norges Bank.
function X = xlags(data,p)
%==========================================================================
% inputs:
% data = raw data
% p = lags
%
% outputs:
% X = lagged dependent variables
%==========================================================================

[n,m] = size(data);

X = zeros(n-p,m*p);

for ii = 1:p
    
    X(:,(1:m)+(ii-1)*m) = data(p-ii+1:end-ii,:);
    
end

