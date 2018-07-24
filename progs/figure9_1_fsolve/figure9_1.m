% FIGURE9_1.M 
%
% Kilian and Lutkepohl (2017), Structural VAR Analysis, Cambridge University Press.
% This file generates Figure 9.1

clear

global A SIGMA p

data; [t,q]=size(y); p=4; h=12;
[A,SIGMA,U,V]=olsvarc(y,p);		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1 of method of moments
% Set seed (to ensure replicability of fsolve results)
randn('seed',1) 

% Set some options for fsolve
warning off
options=optimset('TolX',1e-4,'TolFun',1e-4,'MaxFunEvals',1e+50000);

% Standard identification based on recursive short-run restrictions
OUTPUT=fsolve('version1',[eye(q) ones(q,1)],options); 
B0inv=inv(OUTPUT(1:q,1:q))*diag(sqrt(OUTPUT(1:q,q+1)),0);
%disp([B0inv chol(SIGMA(1:q,1:q))'])
disp(B0inv)
%B0=OUTPUT(1:q,1:q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 2 of method of Moments
% Set seed (to ensure replicability of fsolve results)
randn('seed',1) 

% Set some options for fsolve
warning off
options=optimset('TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e+10);

% Standard identification based on recursive short-run restrictions
B0inv=fsolve('version2',randn(q,q),options); 
%disp([B0inv chol(SIGMA(1:q,1:q))'])

% If necessary, flip the signs as needed, so diagonal elements of B0inv are 
% positive, e.g.:
B0inv(:,2)=-B0inv(:,2);
disp(B0inv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot structural impulse responses
IRF=irfvar(A,B0inv,p,h);
subplot(1,3,1); plot(0:1:h,cumsum(IRF(1,:)),0:1:h,zeros(1,h+1),'linewidth',3); axis([0 h  0 20])
title('Real Price of Oil','fontsize',16)
subplot(1,3,2); plot(0:1:h,cumsum(IRF(2,:)),0:1:h,zeros(1,h+1),'linewidth',3); axis([0 h -1 1])
title('GDP Deflator','fontsize',16)
subplot(1,3,3); plot(0:1:h,cumsum(IRF(3,:)),0:1:h,zeros(1,h+1),'linewidth',3); axis([0 h -1 1])
title('Real GDP','fontsize',16)


