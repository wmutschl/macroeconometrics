function [yc,ytr]=bnDecomp(y,p)
% PURPOSE: Compute Beveridge-Nelson Decomposition of trend and cycle of
% AR(p) model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: 
%
% y = (t x 1) vector of observable series in levels (or log levels)
%
% p = number of lags in AR(p) models
%
% Output:
%
% yc = cycle component of y
%
% ytr = trend component of y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% Initial transformations
%TLe=size(y,1);
yd=diff(y);
%% Do AR(p) estimation
% Construct variables for regression
yl=latMlag(yd,p);
yl=yl(1+p:end,:);
T=size(yd,1);
c=ones(T-p,1);
x=[c yl];
yy=yd(1+p:end);
% OLS
beta=x\yy;    
%% BN 
% Get companion form, i.e., write AR(p) as AR(1)
[betac,~]=varGetCompForm(beta(2:end)',[],p,1); 

c1=betac/(eye(p)-betac); % same as: betacxx (I-betac)^{-1}
ydd=([zeros(1+p,1);yd]-beta(1));
% Construct matrix of historical lags (nlag x TLe), i.e. conditioning 
% forecasting vector for each t
yD=[];
for i=1:p
    yD=cat(1,yD,(ydd(1+i:end-(p)+i)'));
end;
yD=flipud(yD);

% Selection vector
sel=zeros(1,p);sel(1)=1;
% Compute trend and cycle
ytr=y+(sel*c1*yD)';
yc=(sel*c1*yD)';


