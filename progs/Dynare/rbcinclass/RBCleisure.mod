%% Declare Variables and Parameters
var Y C K L A R W I growth;
varexo eps_A;
parameters alph betta delt gam pssi rhoA sigA etaC etaL;
varobs growth;
%% Calibration of parameters
K_o_Y = 10; % avg capital productivity
I_o_Y = 0.25; % avg investment ouput ratio
alph = 0.35;  % cobb douglas
gam  = 1; % normalize
Lbar = 1/3;
rhoA = 0.9;
sigA = 0.6;
etaC = 1;    % risk aversion of households
etaL = 1;  % inverse of intertemporal elasticity of substitution w.r.t leisure
             % Frish elasticity of labor

delt = I_o_Y / K_o_Y;
betta = 1/(alph/K_o_Y+1-delt);
Abar = 1;
Rbar = 1/betta+delt-1;
K_o_L = ((alph*Abar)/Rbar)^(1/(1-alph));
Kbar = K_o_L*Lbar;
Ybar = Kbar/K_o_Y;
Ibar = delt*Kbar;
Wbar = (1-alph)*Abar*K_o_L^alph;
Cbar = Ybar - Ibar;
pssi = gam*(Cbar/Lbar)^(-etaC)*Wbar*(Lbar/(1-Lbar))^(-etaL);

%% Model equations
model;
% Euler equation
gam*C^(-etaC) = betta*gam*C(+1)^(-etaC)*(1-delt+R(+1));
% Labor supply
W = - (-pssi*(1-L)^(-etaL))/(gam*C^(-etaC));
% Capital demand
R = alph*Y/K(-1);
% Labor demand
W= (1-alph)*Y/L;
% Production function
Y = A*K(-1)^alph*L^(1-alph);
% Capital accumulation
K = (1-delt)*K(-1) + I;
% Goods market Clearing
Y = C + I;
% Technology process
log(A) = rhoA*log(A(-1)) + sigA*eps_A;
% Growth measurement equation
growth = (Y-Y(-1))/Y(-1);
end;

%% Shock covariance declaration
shocks;
var eps_A = 1;
end;

%% Steady State (for etaC=etaL=1)
steady_state_model;
A = 1;
R = 1/betta+delt-1;
K_L = ((alph*A)/R)^(1/(1-alph));
W = (1-alph)*A*K_L^alph;
I_L = delt*K_L;
Y_L = A*K_L^alph;
C_L = Y_L - I_L;
L = gam/pssi*C_L^(-1)*W/(1+gam/pssi*C_L^(-1)*W);
C=C_L*L;
I=I_L*L;
K=K_L*L;
Y=Y_L*L;
growth = 0;
end;

%% Computations
resid;
steady(solve_algo=4);
resid;
stoch_simul(order=1,irf=0,periods=10000);
save('simdat','growth');

% bayesian syntax:    
estimated_params;
% parameter namen, startwert, untere Grenze, obere Grenze;
  	alph, normal_pdf, .3, .1, 0.1,0.9;
%   betta, uniform_pdf, 0.98, 0.999, 0.9, 0.9999;
%   delt, uniform_pdf, 0.1, 0.3, 0.05, 0.5;
%  	theta, beta_pdf, .35, .05;
%  	pssi, gamma_pdf, 2,.5, 0.5, 6;
    rhoA, beta_pdf, 0.8,0.1,  0.01, 0.99;
    sigA, inv_gamma_pdf, 0.5, 4, 0.0001,5;
end;

%estimated_params;
%alph, 0.4, 0.0001, 0.9999;
%betta, 0.99, 0.8, 0.99999;
%delt, 0.1, 0.0001, 0.9999;
%gam, 1, 0.5, 3;
%pssi, 1.7, 1, 2;
%rhoA, 0.8, 0.1, 0.9999;
%sigA, 0.6, 0.01, 1.5;
%stderr eps_A, 0.6, 0.01, 1.5;
%end;
    
estimation(datafile='simdat.mat',
           nobs=200,
           first_obs=203,
           mode_compute=4,
           mode_check,
           mh_replic=0);

estimation(datafile='simdat.mat',
           nobs=200,
           first_obs=203,
           mode_compute=0,
           mode_file = 'RBCleisure_mode',
           mode_check,
           mh_replic=5000,
           mh_nblocks = 2,           
           mh_jscale=1.5);
%steady; check;
