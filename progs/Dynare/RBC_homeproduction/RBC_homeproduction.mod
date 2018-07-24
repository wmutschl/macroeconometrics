var Y,Cm,Ch,I, K, Lm,Lh, W,R,A,B;
varexo epsA epsB;
parameters alph, betta, delt, gam, omeg, etta, thet, rhoA, rhoB, sigA, sigB, Abar, Bbar;

alph = 0.35;
betta = 0.97;
delt=0.06;
gam = 0.4;
omeg = 0;
etta=0.8;
thet =1;
rhoA = 0.95;
sigA = 0.01;
rhoB = 0.95;
sigB = 0.01;
Abar = 1;
Bbar = 1;

model(block,bytecode);
betta*Cm(+1)^(etta-1)/(omeg*Cm(+1)^etta +(1-omeg)*Ch(+1)^etta)*(R(+1)+1-delt) = Cm^(etta-1)/(omeg*Cm^etta +(1-omeg)*Ch^etta);
(1-gam)/(1-Lm-Lh) = gam*omeg*Cm^(etta-1)/(omeg*Cm^etta +(1-omeg)*Ch^etta)*W;
(1-gam)/(1-Lm-Lh) = gam*(1-omeg)*thet*Ch^etta/((omeg*Cm^etta +(1-omeg)*Ch^etta)*Lh);
W = (1-alph)*A*(K(-1)/Lm)^alph;
R = alph*A*(Lm/K(-1))^(1-alph);
Ch = B*Lh^thet;
Y = A*K(-1)^alph*Lm^(1-alph);
K = (1-delt)*K(-1) + I;
Y = Cm+I;
log(A) = (1-rhoA)*log(Abar)+rhoA*log(A(-1)) + epsA;
log(B) = (1-rhoB)*log(Bbar)+rhoB*log(B(-1)) + epsB;
end;

initval;
epsA=0;
epsB=0;
A=1;
B=1;
Y=0.459779;
Cm=0.353592;
Ch=0.199197;
Lm=0.22251;
Lh=0.133077;
K=1.76978;
I=0.106187;
W=(1-alph)*Y/Lm;
R=alph*Y/K;
end;
steady;

%homotopy_setup;
%omeg, 0.5,0.9999;
%end;
%steady(homotopy_mode=1,homotopy_steps = 5,solve_algo = 5, markowitz = 5);

% stoch_simul(order=1);
%     save_params_and_steady_state('ss1.txt');

