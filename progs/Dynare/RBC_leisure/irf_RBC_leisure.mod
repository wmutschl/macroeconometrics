% =========================================================================
% Impulse responses for
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in either CRRA or log utility function, 
%   dependent on Dynare macro CRRA (see below)
% - Cobb Douglas production function with labor & capital in a perfect
%   competition setting
% - Total Factor Productivity shock
% -------------------------------------------------------------------------
% Based on codes by Martin Andreasen and Johannes Pfeifer
% Willi Mutschler, February 2018
% willi@mutschler.eu
% =========================================================================

% =========================================================================
% Dynare preprocessor
% =========================================================================
% Set options for Dynare preprocessor:
%  - CRRA: specification of utility function: 1 for CRRA, else for log
%  - elaboratecalib: 0 for standard values or 1 for elaborate way to 
%                    calibrate parameters based on steady state relations
%  - irfhorizon: horizon of impulse response function
% =========================================================================
@#define CRRA = 0
@#define elaboratecalib = 1
@#define irfhorizon = 40

% =========================================================================
% Declare variables, shocks and parameters
% =========================================================================
@#include "RBC_leisure_symdecls.mod"

% =========================================================================
% Calibrate parameter values based on OECD data
% =========================================================================
@# if elaboratecalib == 1
    @#include "RBC_leisure_params_elaborate.mod"
@# else
    @#include "RBC_leisure_params_simple.mod"
@# endif

% =========================================================================
% Model equations
% =========================================================================
@#include "RBC_leisure_modeqs.mod"

% =========================================================================
% Steady state computations either initval or steady_state_model
% dependent on form of utility function
% =========================================================================
@# if CRRA 
    @#include "RBC_leisure_initval.mod"
@# else
    @#include "RBC_leisure_ssmodel.mod"
@# endif

% =========================================================================
% Declare settings for shocks
% =========================================================================
@#include "RBC_leisure_epsA.mod"

% =========================================================================
% Computations
% =========================================================================
steady; % compute steady state given the starting values
resid;  % check the starting values for the steady state
check;  % check Blanchard-Khan conditions

% Approximate policy functions with first-order perturbation methods
% Compute and plot impulse response function set above in irfhorizon
% The irf functions are stored in variablename_shockname or in oo_.irfs
stoch_simul(order=1,irf=@{irfhorizon}); 
