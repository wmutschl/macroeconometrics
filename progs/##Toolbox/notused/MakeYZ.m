% function [Z,Y] = MakeZY(p,y,ydettrend)
% K = size(y,2);
% T = size(y,1);
% init = 0;
% Z = ones(K*p,T)*init;
% for ii=1:p
%     Z((K*(ii-1)+1):K*ii,1+ii:T)=y(1:T-ii,:)';
% end
% Z = [ydettrend; Z];        
% 
% Z=Z(:,1+p:end);
% Y=y(1+p:end,:)';


function [Y, Z] = MakeYZ(DATA,lags,const)
% =======================================================================
% Builds the VAR process from the data-matrix DATA. It orders the data into
% the Y and X matrix --> Example: [x y] = [x(-1) y(-1) x(-2) y(-2)]
% =======================================================================
% [Y, X] = VARmakexy(DATA, lags, const)
% -----------------------------------------------------------------------
% INPUT
%   DATA: matrix containing the original data
%   lags: lag order of the VAR
%   const : 0, no constant, no trend
%           1, constant, no trend
%           2, constant, trend
% -----------------------------------------------------------------------
% OUTPUT
%   Y: dependent variable
%   Z: independent variable
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


% Get dimesion of DATA
nobs = size(DATA,2);

% Y matrix 
Y = DATA(:,lags+1:end);

% Z-matrix 
Z=[];
for jj=0:lags-1
    Z = [DATA(:,jj+1:nobs-lags+jj); Z];
end

if const==1 %constant        
       Z = [ones(1,nobs-lags); Z];       
elseif const==2 % time trend and constant
        trend=1:size(Z,2);
        Z = [ones(1,nobs-lags); trend; Z];     
end
