var Y,C,Ci,Cj,I,Ii, K,Ki, L,Li,Lj, W,R,A;
varexo epsA;
parameters alph, betta, delt, gam, omeg, rhoA, sigA;

alph = 0.35;
betta = 0.97;
delt=0.06;
gam = 0.4;
omeg = 1;
rhoA = 0.95;
sigA = 0.01;
phhi = 0;

model;
Ci(+1)/Ci = betta*(1-delt+R(+1));
Cj = W*Lj;
Ci = (gam/(1-gam))*(1-Li)*W;
Cj = (gam/(1-gam))*(1-Lj)*W;
L = omeg*Li + (1-omeg)*Lj;
C = omeg*Ci + (1-omeg)*Cj;
K = omeg*Ki;
I = omeg*Ii;

W = (1-alph)*A*(K(-1)/L)^alph;
R = alph*A*(L/K(-1))^(1-alph);
Y = A*K(-1)^alph*L^(1-alph);
Ki = (1-delt)*Ki(-1) + Ii;
Y = C+I;
log(A) = rhoA*log(A(-1)) + epsA;
end;

initval;
Y=1;
C=0.8;
Ci=0.6;
Cj=0.2;
L=0.3;
Li=0.3;
Lj=0.3;
K=3.5;
I=0.2;
Ii = 0.3;
W=(1-alph)*Y/L;
R=alph*Y/K;
A=1;
epsA=0;
end;

steady;
check;

shocks;
var epsA = sigA^2;
end;

stoch_simul(order=1) Y C Ci Cj I K L Li Lj W;