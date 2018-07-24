% =========================================================================
% Calibrate parameter values based on steady state relationships for:
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in either CRRA or log utility function, 
%   dependent on Dynare macro CRRA.
% - Cobb Douglas production function with labor & capital in a perfect
%   competition setting
% - Total Factor Productivity shock
% -------------------------------------------------------------------------
% Based on codes by Martin Andreasen and Johannes Pfeifer
% Willi Mutschler, February 2018
% willi@mutschler.eu
% =========================================================================
K_o_Y = 10;   % average capital capital productivity K/Y
I_o_Y = 0.25; % average investment output ratio I/Y
alph  = 0.35; % capital share based
gam   = 1;    % normalize consumption preference parameter to unity
Lbar  = 1/3;  % steady state labor share
rhoA  = 0.9;  % technology autocorrelation based on linearly detrended Solow residual
sigA  = 0.6;  % technology standard deviation  based on linearly detrended Solow residual
delt  = I_o_Y/K_o_Y;             % quarterly depreciation rate
betta = 1/(alph/K_o_Y+1-delt); % discount factor
@# if CRRA
    etaC  = 2;    % risk aversion of households
    etaL  = 1.5;  % inverse of intertemporal elasticity of substitution w.r.t leisure
@#endif
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
@#if CRRA
    pssi = gam*(Cbar/Lbar)^(-etaC)*Wbar*(Lbar/(1-Lbar))^(-etaL);
@#else
    pssi = gam*(Cbar/Lbar)^(-1)*Wbar*(Lbar/(1-Lbar))^(-1);    
@#endif
