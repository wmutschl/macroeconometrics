var Y,C,I,K,L,W,R,A;
varexo epsA;
parameters alph, betta, delt, gam, rhoA, phhi,sigA;

alph = 0.35;
betta = 0.97;
delt=0.06;
gam = 0.4;
rhoA = 0.95;
sigA = 0.01;
phhi = 0;

model;
#cback = 1/(C -phhi*C(-1));
#c = 1/(C(+1) -phhi*C);
#cp = 1/(C(+2) - phhi*C(+1));
(cback-betta*phhi*c)/(c-betta*phhi*cp) = betta*(R(+1)+1-delt);
(gam*cback - betta*gam*phhi*c)*W = (1-gam)/(1-L);

% (1/(C(-1) -phhi*C(-2))-betta*phhi*1/(C -phhi*C(-1)))/(1/(C -phhi*C(-1))-betta*phhi*1/(C(+1) - phhi*C)) = betta*(R(+1)+1-delt);
% (gam*1/(C -phhi*C(-1)) - betta*gam*phhi*1/(C(+1) - phhi*C))*W = (1-gam)/(1-L);
W = (1-alph)*A*(K(-1)/L)^alph;
R = alph*A*(L/K(-1))^(1-alph);
Y = A*K(-1)^alph*L^(1-alph);
K = (1-delt)*K(-1) + I;
Y = C+I;
log(A) = rhoA*log(A(-1)) + epsA;
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
epsA=0;
end;

steady;
check;

shocks;
var epsA = sigA^2;
end;

stoch_simul(order=1);