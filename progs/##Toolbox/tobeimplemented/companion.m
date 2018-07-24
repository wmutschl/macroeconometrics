%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                %
%          MatLab function for estimating eigenvalues            %
%          of the companion matrix                               %
%                                                                %
%          Function: companion                                   %
%                                                                %
%          Input:                                                %
%          	            Z0: nxT matrix                           %
%                       Z1: nq1xT matrix                         %
%                       Z2: n(p-1)+d matrix                      %
%                       r: cointegration rank                    %
%                       p: lag length                            %
%                       n: number of endogenous variables        %
%                       q1: number of exogenous I(1) variables   %
%                       q0: number of exogenous I(0) variables   %
%                       d: number of deterministic variables     %
%                       T: sample size                           %
%                                                                %
%          Output:                                               %
%          	            eigenvalues: nx1 vector of eigenvalues   %
%                           of the companion matrix              %
%                                                                %
%          This version:        January 17, 2005                 %
%          Written by: 	        Anders Warne                     %
%                               Copyright (c) 2001-2005          %
%                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigenvalues = companion(Z0,Z1,Z2,r,p,n,q1,q0,d,T)

[S_00,S_01,S_10,S_11,M_02,M_22,M_12,beta,lambda,trace_test] = mlcointbasic(Z0,Z1,Z2,r);

alpha = (S_01*beta)/(beta'*S_11*beta);
if isempty(Z2)==0;
   Psi = (M_02/M_22)-(alpha*beta'*(M_12/M_22));
else;
   Psi = [];
end;

Pi = alpha*beta(1:n,:)';
if p>=2;
   Phi = Psi(:,d+1:d+n*(p-1));
end;

ip = zeros(p,p);
i = 1;
while i<=p;
   j = 1;
   while j<=p;
      if j>=i;
         ip(j,i)=1;
      end;
      j = j+1;
   end;
   i = i+1;
end;

%
% Computing levels VAR coefficients from
% VECM coefficients
%

ip = kron(ip,eye(n));
if p>=2;
   M = [(Pi+eye(n)) -Phi]/ip;
   M = [M;[eye(n*(p-1)) zeros(n*(p-1),n)]];
else;
   M = (Pi+eye(n))/ip;
end;

eigenvalues = eig(M);

%
% end of function companion 
%
