% =========================================================================
% Calibrate parameter values given standard values for:
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
alph  = 0.35;
betta = 0.99;
delt  = 0.025;
gam   = 1;
pssi  = 1.6;
rhoA  = 0.9;
sigA  = 0.6;
@# if CRRA
    etaC  = 2;    % risk aversion of households
    etaL  = 1.5;  % inverse of intertemporal elasticity of substitution w.r.t leisure
@#endif
