% =========================================================================
% Declare enodgenous and exogenous variables as well as parameters for:
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
var 
    Y        ${Y}$  (long_name='output')
    C        ${C}$  (long_name='consumption')
    K        ${K}$  (long_name='capital')
    L        ${L}$  (long_name='hours')
    A        ${Z}$  (long_name='total factor productivity')
    R        ${R}$  (long_name='interest rate')
    W        ${W}$  (long_name='real wage')
    I        ${I}$  (long_name='investment')
;

varexo
    eps_A    ${\varepsilon_A}$ (long_name='TFP shock')
;

parameters
    alph    ${\alpha}$   (long_name='capital share')
    betta   ${\beta}$    (long_name='discount factor')
    delt    ${\delta}$   (long_name='depreciation rate')
    gam     ${\gamma}$   (long_name='consumption utility parameter')
    pssi    ${\psi}$     (long_name='labor disutility parameter')
    rhoA    ${\rho_A}$   (long_name='persistence TFP')
    sigA    ${\sigma_A}$ (long_name='standard deviation TFP shock')
    @#if CRRA
    etaC    ${\eta_C}$   (long_name='risk aversion of household')
    etaL    ${\eta_L}$   (long_name='inverse of labor elasticity of substitution')
    @#endif
;