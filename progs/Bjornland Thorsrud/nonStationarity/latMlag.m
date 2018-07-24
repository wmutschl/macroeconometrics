function xlag=latMlag(x,p,init)
% PURPOSE: Generate lags of matrix (vector), following common VAR
% standards, eg: x1(t-1) x2(t-1) ... xN(t-1)... x1(t-p)... xN(t-p). 
% NOTE: This is NOT as done in the Lesage package, which uses
% x1(t-1)...x1(t-p) ... xN(t-1)...xN(t-1)..xN(t-p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
%   x = matrix or vector, (t x n) where t is number of obs, and n is number
%   of variables.
%
%   p = integer. Number of lags. 
%
%   init = (optional) scalar value to feed initial missing values
%   (default = 0)
%
% Output:
%
%   xlag = matrix or vector with lags of x. The output has the same size as
%   the input matrix or vector. Missing values at the beginning (due to
%   taking lags) have been replaced by the init argument. 
%
% Usage:
%
%   xlag=latMlag(x,p,init)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% leifath@gmail.com
% 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ==1 
    p = 1; % default value
    init = 0;
elseif nargin == 2
    init = 0;
end;

if nargin > 3
    error('mlag:err','Wrong # of input arguments');
end;

[nobs, nvar] = size(x);

xlag = ones(nobs,nvar*p)*init;

for ii=1:p
    xlag(1+ii:nobs,(nvar*(ii-1)+1):nvar*ii)=x(1:nobs-ii,:);
end
