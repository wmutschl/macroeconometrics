% =========================================================================
% Compute steady state using a steady_state_model block for
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in log utility function (not CRRA)
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
% Note that I use DYNARE's special macro language (@#)
% That is, whenever you see " @# " the DYNARE preprocessor will 
% include source files, replicate blocks of equations through loops (for), 
% conditionally execute some code (if/then/else) or substitute expressions
% -------------------------------------------------------------------------
% Set options for Dynare preprocessor:
%  - CRRA: specification of utility function: 1 for CRRA, else for log
%  - elaboratecalib: 0 for standard values or 1 for elaborate way to 
%                    calibrate parameters based on steady state relations
% =========================================================================
@#define CRRA = 0
@#define elaboratecalib = 1

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
% Steady state using a steady_state_model block
% =========================================================================
@#include "RBC_leisure_ssmodel.mod"

% =========================================================================
% Computations
% =========================================================================
steady; % compute steady state given the starting values
resid;  % check the starting values for the steady state
