function mse=getMse(betac,sigma,h)
% PURPOSE: Calculate MSE from VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% betac = coefficient matrix. Companion form
%
% sigma = Covariance matrix of VAR residuals, (n x n) 
%
% h = Scalar. Number of horizons 
%
% Output:
%
% mse = mse array. (n x n x h)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Placeholders and empty output
n=size(sigma,1);
nb=size(betac,1);
mseb=zeros(nb,nb);
mse=zeros(n,n,h);

% Selection matrix
sel=[eye(n) zeros(n,nb-n)];

% Put sigma into the upper left corner of the bigger matrix
mseb(1:n,1:n)=sigma;
mseb0=mseb;
% First observation is just equal to sigma (B_{1}=I)
mse(:,:,1)=sigma;
for i=2:h
    b=betac^(i-1);    
    mseb=mseb+b*mseb0*b';
    % Only extract the current period
    mse(:,:,i)=sel*mseb*sel';
end;






