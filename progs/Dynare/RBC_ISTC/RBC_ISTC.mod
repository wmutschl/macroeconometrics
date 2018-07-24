var Y,C,I,K,L,W,R,A,Z;
varexo epsA epsZ;
parameters alph, betta, delt, gam, rhoA, rhoZ, sigA, sigZ;

alph = 0.35;
betta = 0.97;
delt=0.06;
gam = 0.4;
rhoA = 0.95;
sigA = 0.01;
rhoZ = 0.95;
sigZ = 0.01;

model;
C(+1)/C = betta*Z/Z(+1)*(1-delt+R(+1));
C = (gam/(1-gam))*(1-L)*W;

W = (1-alph)*A*(K(-1)/L)^alph;
R = alph*A*(L/K(-1))^(1-alph);
Y = A*K(-1)^alph*L^(1-alph);
K = (1-delt)*K(-1) + I;
Y = C+I;
log(A) = rhoA*log(A(-1)) + epsA;
log(Z) = rhoZ*log(Z(-1)) + epsZ;
end;

initval;
Y=1;
C=0.8;
L=0.3;
K=3.5;
I=0.2;
W=(1-alph)*Y/L;
R=alph*Y/K;
A=1;
Z=1;
epsA=0;
epsZ=0;
end;

steady;
check;

shocks;
var epsA = sigA^2;
var epsZ = sigZ^2;
end;

stoch_simul(order=1);