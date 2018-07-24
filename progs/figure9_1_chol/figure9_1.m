% FIGURE9_1.M 
%
% Kilian and Lutkepohl (2017), Structural VAR Analysis, Cambridge University Press.
% This file generates Figure 9.1

clear

global h q 

data; [t,q]=size(y); p=4; h=12;
[A,SIGMA,U,V]=olsvarc(y,p);		
B0inv=chol(SIGMA(1:q,1:q))';

h=12; IRF=irfvar(A,B0inv,p,h);
subplot(1,3,1); plot(0:1:h,cumsum(IRF(1,:)),0:1:h,zeros(1,h+1),'k-','linewidth',3); axis([0 h  0 20])
title('Real Price of Oil','fontsize',16)
xlabel('Quarters')
subplot(1,3,2); plot(0:1:h,cumsum(IRF(2,:)),0:1:h,zeros(1,h+1),'k-','linewidth',3); axis([0 h -1 1])
title('GDP Deflator','fontsize',16)
xlabel('Quarters')
subplot(1,3,3); plot(0:1:h,cumsum(IRF(3,:)),0:1:h,zeros(1,h+1),'k-','linewidth',3); axis([0 h -1 1])
title('Real GDP','fontsize',16)
xlabel('Quarters')

