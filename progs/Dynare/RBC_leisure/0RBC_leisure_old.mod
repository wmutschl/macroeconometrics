% =========================================================================
% Baseline RBC model with TFP shock, consumption and leisure in CRRA
% utility function and Cobb Douglas production function with labor & capital
% Computations:
% c) (log as special case)
% a) 3 ways to compute steady-state:
%   - initival
%   - external function RBC_leisuresteadystate.m
%   - steady_state_model (only with log utility)
% b) 2 ways to calibrate parameters
%    - elaborate way using steady state relationships
%    - standard calibration used in literature
% Based on codes by Martin Andreasen and Johannes Pfeifer
% =========================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =========================================================================

% =========================================================================
%   Note that I use DYNARE special language (@#), so we do not need to uncomment stuff below
%   That is, whenever you see " @# " the DYNARE preprocessor will execute only
%   those computations depending on the macro variable scenario
% Set options for Dynare preprocessor:
%  - logutility: specification of utility function: 1 for log, else CRRA
%  - steadystatemethod: 1 initval, 2 steady_state_model (only for log utility)
%                       or else using external matlab file
%    Note that for initval you need to rename the external matlab file
%    RBC_leisure_steadystate.m to something else
%  - elaboratecalib: calibration of parameters: 1 for elaborate way, else
%                    standard values
% Set macro directive for which scenario (this is just for convenience), no semicolon
% =========================================================================
@#define logutility = 1
@#define steadystatemethod = 2
@#define elaboratecalib = 1
% scenario = 1: Compute steady state
% scenario = 2: Stochastic model - temporary shock on technology
% scenario = 3: Stochastic model - simulate data
% scenario = 4: Maximum likelihood estimation given a large sample size
% scenario = 5: Maximum likelihood estimation given a small sample size
% scenario = 6: Bayesian Estimation given a large sample size
% scenario = 7: Bayesian Estimation given a small sample size
@#define scenario = 4

% =========================================================================
% Declare variables, shocks and parameters
% =========================================================================
var 
    Y        ${Y}$  (long_name='output')
    C        ${C}$  (long_name='consumption')
    K        ${K}$  (long_name='capital')
    L        ${L}$  (long_name='hours')
    A        ${Z}$  (long_name='total factor productivity')
    R        ${R}$  (long_name='interest rate')
    W        ${W}$  (long_name='real wage')
    I        ${I}$  (long_name='investment') 
    dy       ${dy}$ (long_name='output growth')
;

varobs dy L;

varexo
    eps_A    ${\varepsilon_A}$ (long_name='TFP shock')
;

    
parameters
    alph    ${\alpha}$   (long_name='capital share')
    betta   ${\beta}$    (long_name='discount factor')
    delt    ${\delta}$   (long_name='depreciation rate')
    etaC    ${\eta_C}$   (long_name='risk aversion of household')
    etaL    ${\eta_L}$   (long_name='inverse of labor elasticity of substitution')
    gam     ${\gamma}$   (long_name='consumption utility parameter')
    pssi    ${\psi}$     (long_name='labor disutility parameter')
    rhoA    ${\rho_A}$   (long_name='persistence TFP')
    sigA    ${\sigma_A}$ (long_name='standard deviation TFP shock')
;

% =========================================================================
% Calibrate parameter values based on OECD data
% =========================================================================
@# if logutility == 1
    etaC  = 1;    % risk aversion of households
    etaL  = 1;    % inverse of intertemporal elasticity of substitution w.r.t leisure
@# else
    etaC  = 2;    % risk aversion of households
    etaL  = 1.5;  % inverse of intertemporal elasticity of substitution w.r.t leisure
@#endif

@# if elaboratecalib == 1
    K_o_Y = 10;   % average capital capital productivity K/Y
    I_o_Y = 0.25; % average investment output ratio I/Y
    alph  = 0.35; % capital share based
    gam   = 1;    % normalize consumption preference parameter to unity
    Lbar  = 1/3;  % steady state labor share
    rhoA  = 0.9;  % technology autocorrelation based on linearly detrended Solow residual
    sigA  = 0.6;  % technology standard deviation  based on linearly detrended Solow residual
    delt  = I_o_Y/K_o_Y;             % quarterly depreciation rate
    betta = 1/(alph/K_o_Y+1-delt); % discount factor
    Abar = 1;                                  % steady state A
    Rbar = 1/betta+delt-1;                     % steady state R
    K_o_L = (alph*Abar/Rbar)^(1/(1-alph));     % steady state K/L
    Kbar = K_o_L*Lbar;                         % steady state K
    Ybar = Kbar/K_o_Y;                         % steady state Y
    Ibar = delt*Kbar;                          % steady state I
    Ybar = Abar*Kbar^alph*Lbar^(1-alph);       % steady state Y
    Wbar = (1-alph)*Abar*(K_o_L)^alph;         % steady state W
    Cbar = Ybar - Ibar;                        % steady state C
    % labor preference parameter given from steady state relationship
    pssi = gam*(Cbar/Lbar)^(-etaC)*Wbar*(Lbar/(1-Lbar))^(-etaL);
@# else
    alph  = 0.35;
    betta = 0.99;
    delt  = 0.025;
    gam   = 1;
    pssi  = 1.6;
    rhoA  = 0.9;
    sigA  = 0.6;
@# endif
sigmeas = 0.1;
% =========================================================================
% Model equations
% =========================================================================
model;
% auxiliary expressions, these will be substituted by Dynare preprocessor
#UC  = gam*C^(-etaC);
#UCp = gam*C(+1)^(-etaC);
#UL  = -pssi*(1-L)^(-etaL);

[name='Euler equation']
UC=betta*UCp*(1-delt+R);

[name='labor supply']
W=-UL/UC;

[name='capital law of motion'] 
K=(1-delt)*K(-1)+I;

[name='market clearing/resource constraint']
Y=I+C;

[name='production function']
Y=A*K(-1)^alph*L^(1-alph);

[name='labor demand']
W=(1-alph)*Y/L;

[name='capital demand']
R=alph*Y/K(-1);

[name='exogenous TFP process']
log(A)=rhoA*log(A(-1))+eps_A;

[name='output growth']
dy = Y/Y(-1)-1 + meas;
end;

% =========================================================================
% Define shock covariance matrix
% =========================================================================
shocks;
    var eps_A = sigA^2;
end;

% =========================================================================
% Define steady state computation
% =========================================================================
% Analytical steady state using steady_state_model for log-utility
@# if steadystatemethod == 2 && logutility == 1
steady_state_model;
    A = 1;                             % technology innovations
    R = 1/betta+delt-1;                % rental rate of capital
    K_L = ((alph*A)/R)^(1/(1-alph));   % capital divided by labor
    W = (1-alph)*A*K_L^alph;           % wage level
    I_L = delt*K_L;                    % investment divided by labor
    Y_L = A*K_L^alph;                  % output divided by labor
    C_L = Y_L-I_L;                     % consumption divided by labor
    % closed-form expression for labor when using log utility
    L = gam/pssi*C_L^(-1)*W/(1+gam/pssi*C_L^(-1)*W);
    C = C_L*L; % consumption level
    I = I_L*L; % investment level
    K = K_L*L; % capital level
    Y = Y_L*L; % output level
    dy = 0;    % output growth
end;

@# else
    @# if steadystatemethod == 1
% Numerical steady state: provide initial values for optimizer
% Note that you need to rename RBC_leisure_steadystate.m otherwise it will
% be used by Dynare instead
initval;
    A =  1;
    Y =  1;
    C =  0.8;
    K =  10;
    L =  1/3;
    R =  0.03;
    W =  2.25;
    I =  0.28;
    dy = 0;
end;
resid;  % check equations for initial values if optimizer fails

    @# else
% Steady state is computed by external matlab file RBC_leisure_steadystate.m
    @#endif
@#endif

% =========================================================================
% SCENARIO 1: Compute steady state and check Blanchard-Khan conditions
% =========================================================================
@#if scenario >= 1
    steady; % compute steady state given the starting values
    resid;  % check the starting values for the steady state
    check;  % check Blanchard-Khan conditions
@#endif

% =========================================================================
% SCENARIO 2: Stochastic model - temporary shock
% =========================================================================
@#if scenario == 2    
    stoch_simul(order=1,nocorr,nomoments,irf=40);
@#endif

% =========================================================================
% SCENARIO 3: Stochastic model - simulate data
% =========================================================================
@#if scenario == 3
    % Simulate 1000 periods
    stoch_simul(order=1,nocorr,nomoments,irf=0,periods=10000);
    % Save data to mat file
    save('simdat','dy',');
@#endif


% =========================================================================
% SCENARIO 4: Maximum likelihood estimation given a large sample size
% =========================================================================
@#if scenario == 4
    % maximum likelihood estimation syntax:    
    estimated_params;
    % parameter, initial value, lower bound, upper bound;
    alph,        0.33,          0.01,        0.99;
    %betta,       0.99,          0.8,         0.999999;
    %delt,        0.02,          0.001,       0.5;
    gam,  1, 0.1, 2;
    %pssi,  1.5, 0.1,   5;
    rhoA,       0.9,            0.5,         0.9999999;
    stderr eps_A,       0.5,            0.01,        1; 
    @#if logutility == 0
    etaC, 
    %etaL,
    @#endif
    end;
    estimation(order=1,datafile='simdat.mat',nobs=10000,first_obs=1,mode_compute=4,mode_check);
@#endif




% =========================================================================
% SCENARIO 5: Maximum Likelihood Estimation small sample size
% =========================================================================
@#if scenario == 5
    shocks;
        var eps_a; stderr 1;
        var eps_tau; stderr 1;     
    end;
    % Maximum-Likelihood-Estimation:
    % parameter, lower bound, upper bound
    estimated_params;
        rho,      0.9, 0.001, 0.99;
        lambda,   0.5, 0.001, 0.99;  
    end;
    estimation(order=1,datafile=simdat,nobs=100,first_obs=51,nograph);        
@#endif




% =========================================================================
% SCENARIO 6: Bayesian Estimation small sample size
% =========================================================================
@#if scenario == 6
   shocks;
        var eps_a; stderr 1;
        var eps_tau; stderr 1;     
    end;
    %Bayesian estimation
    %Syntax: parameter, prior type, prior par1, prior par2,... 
    estimated_params;
       rho, beta_pdf, 0.9, .04;
       lambda, normal_pdf, 0.5, .04; 
    end;
    estimation(datafile=simdat,first_obs=51,nobs=100,mode_check,mh_replic=15000,mh_jscale=2.5,mh_nblocks=2);
@#endif
