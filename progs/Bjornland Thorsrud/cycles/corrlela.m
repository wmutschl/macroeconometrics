function [x,p]=corrlela(X,maxLeadLag)
% PURPOSE: Produce correlations between two variables on both lead and lags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   
%   X = matrix. (t x n), where t is the number of observations and n=2, eg.
%   two variables.
%
%   varargin = optional. If provided string followed by value. One or more
%   of the following
%
%       'maxLeadLag',2 = defines the maximum lead and lag correlations to
%       consider. Default=2. 
%
% Output:
%
%   x = vector. (m x 1), where m is 2*maxLeadLag+1 for the contemperanous
%   correlation. The ordering is as follows ; 
%   [corr(x(1:end-m+1,1),x(m:end,2) x(1:end-m+2,1),x(m+1:end,2) ...etc. ]' 
%
%   p = vector. (m x 1), where size and ordering as in x. Containg p-values
%   for correlations. Values below 0.05 are significant. 
%
% Note: x1=X(:,1) x2=X(:,2), output x if maxLeadLag=2 will be a vector with
% 5 elements: x(1)=corr(x1(1:end-i),x2(1+i-1:end) i.e. x1 leads x2
% x(5)=corr(x1(1+i-1:end),x2(1:end-i) i.e. x1 lags x2
% x(3)= contemperanous correlation
%
%  x2    -------------
%  x1        -------------        This is the correlation coeff at (x(end,1))
%
%  x2        -------------
%  x1    -------------            This is the correlation coeff at (x(1,1))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% 2009-2010
% leifath@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default
options.maxLeadLag=maxLeadLag;

% separate the input 
x1=X(:,1);
x2=X(:,2);

% empty output matrix
x=nan(options.maxLeadLag*2+1,1);
p=x;
% make index to get the correct lead and lag structure
index=(options.maxLeadLag:-1:-(options.maxLeadLag));

for i=1:options.maxLeadLag*2+1
    
    if index(i)==0
        [xi pi]=corr([x1(:,1) x2(:,1)]);
    else
        % x is the correlation coeff, while p is the p-value. A value less than
        % 0.05 makes the correlation significant different from 0.
        if index(i)>0
            [xi pi]=corr([x1(1:end-index(i),1) x2(index(i)+1:end,1)]);
        elseif index(i)<0
            [xi pi]=corr([x1(1-index(i):end,1) x2(1:end+index(i),1)]);
        end;        
    end;
    
    x(i)=xi(1,2);p(i)=pi(1,2);
    
end;       