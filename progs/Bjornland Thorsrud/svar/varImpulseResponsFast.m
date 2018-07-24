function irf=varImpulseResponsFast(betaC,nvar,nlag,hor,A)
% PURPOSE: Faster version (in terms of execution speed) of 
% varImpulseRespons. Useful for e.g. simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   
%   betaC = matrix with the estimated coefficents from the VAR
%   estimation in companion form, e.g: 
%   [b11 b12 ... b1np
%    .     .    .
%    bn1 bn2 ... bnp
%    Inxn      0   ], where n is the number of variables, p is the number
%    of lags and I is the indentity mat. 
%
%   nvar = number of variables, eg. VAR equations. 
%
%   nlag = number of lags used in VAR
%
%   hor = number of response horizons 
%
%   A = matrix, (n x n), where n is number of 
%   variabels. This will typically be the structural contemporaneus matrix 
%   from the SVAR. 
%
% Output:
%   irf = impulse response matrix, (n x n x hor), where n is the number of
%   variables. The first row refers to eqaution 1, and the different
%   columns refers to shocks to variable j. See Hamilton 1994. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% 2010
% leifath@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty psi matrix and irf mat
psi=nan(nvar,nvar,hor); 
irf=psi;

y=zeros(nvar*nlag,nvar*nlag,hor+1);        
y(1:nvar,1:nvar,1)=eye(nvar);
for s=2:hor+1                
    y(:,:,s)=betaC*y(:,:,s-1);
    psi(:,:,s-1)=y(1:nvar,1:nvar,s);    
    irf(:,:,s-1)=psi(:,:,s-1)*A;
end;    
% 15.01.2010: added the contemporaneus (or impact) shocks to the irf. Eg. they start in
% period 0 (not as 1 before) 
irf=cat(3,A,irf(:,:,1:end-1));
