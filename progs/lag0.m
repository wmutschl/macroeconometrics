function out=lag0(x,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For an input vector or matrix x=(a_{1}; a_{2}; ...; a_{n}) where a_{i}%
%are row vectors, it returns the vector or matrix:                     %
%                                                                      % 
%xl0= (0; ... ; 0 ; a_{1}; a_{2}; ....; a_{n-p})                       %
%                                                                      % 
%In other words it lags the variable p periods and places zeros in the % 
%rows corresponding to the first p periods                             % 
%                                                                      % 
%Version is 1.01                                                       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute the number of rows and columns of input x and save in 
%variables R and C respectively
[R,C]=size(x);

%Take the first R-p rows of matrix x
x1=x(1:(R-p),:);
%Preceed them with p rows of zeros and return
out=[zeros(p,C); x1];


%Revision History:
%v 1.01
%Replaced x1=x(1:(length(x)-p),:) with x1=x(1:(R-p),:); 
%Now should work correctly for matrices;