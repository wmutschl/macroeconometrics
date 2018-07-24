function Atilda=longRunIdentification(betac,sigma)
% PURPOSE: Compute simple long run restrictions for structural VAR
% identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
% betac = companion form of VAR coefficient matrix 
%
% sigma = covariance matrix of reduced form residuals (n x n), where n is
% the number of variables in the VAR
%
% Output:
%
% Atilda = structural impact matrix (n x n)
%
% Usage:
%
% Atilda=longRunIdentification(betac,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(sigma,1);

A0=chol(sigma,'lower');    
AA=(eye(size(betac,1))-betac)\eye(size(betac,1));
AA=AA(1:N,1:N);

Cl=AA*A0; 
[q1,r1]=qr(Cl');    
% normalize the diagonal
for ii=1:N;
    if r1(ii,ii)<0
        q1(:,ii)=-q1(:,ii);
    end;
end;                        

Atilda=A0*q1;
    
% alternatice computation
%pp=AA*sigma*AA';                            
%longRunB=chol(pp,'lower');
%Atilda=AA\longRunB;                  