% =========================================================================
% Simulate data for output growth and estimate with Bayesian methods for
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in either CRRA or log utility function, 
%   dependent on Dynare macro CRRA (see below)
% - Cobb Douglas prodduction function with labor & capital in a perfect
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
%  - nperiods: sample size for simulated data
% =========================================================================
@#define CRRA = 0
@#define elaboratecalib = 1
@#define nperiods = 200

% =========================================================================
% Declare variables, shocks and parameters
% Note that we add output growth as additional and observable variable
% =========================================================================
@#include "RBC_leisure_symdecls.mod"
var dy       ${dy}$ (long_name='output growth');
varobs dy;
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
% Note that we add measurement equation for output growth
% =========================================================================
@#include "RBC_leisure_modeqs.mod"
model;
[name='output growth']
dy = Y/Y(-1)-1;
end;
% =========================================================================
% Steady state computations either initval or steady_state_model
% dependent on form of utility function
% Note that we need to add expressions for output growth
% =========================================================================
@# if CRRA 
    @#include "RBC_leisure_initval.mod"
    initval;
    dy=0;
    end;
@# else
    @#include "RBC_leisure_ssmodel.mod"
    steady_state_model;
    dy = 0;
    end;
@# endif

% =========================================================================
% Declare settings for shocks
% =========================================================================
@#include "RBC_leisure_epsA.mod"

% =========================================================================
% Estimation settings for parameter, syntax for Bayesian estimation is:
% parameter name, prior type, prior par1, prior par2,... 
% =========================================================================
estimated_params;
  	alph, normal_pdf, .3, .1, 0.1,0.9;
    rhoA, beta_pdf, 0.8,0.1,  0.01, 0.99;
    stderr eps_A, inv_gamma_pdf, 0.5, 4, 0.0001,5;
end;

% =========================================================================
% Computations
% =========================================================================
steady; % compute steady state given the starting values
resid;  % check the starting values for the steady state
check;  % check Blanchard-Khan conditions

% Approximate policy functions with first-order perturbation methods
% Simulate data for nperiods periods set above
% The simulated data is stored in variablename or in oo_.endo_simul
stoch_simul(order=1,irf=0,periods=@{nperiods}); 
save('simdat','dy'); % Save data to mat file

% Estimate mode of posterior without Metropolis-Hastings algorithm (mh_replic=0)
estimation(datafile='simdat.mat',
           first_obs=1,           % index of first period to use
           nobs=@{nperiods},      % sample size        
           mode_compute=4,        % Choose optimizer (good choices can be 2,4,9,10)
           mode_check,            % plot the posterior for values around the mode for each estimated parameter in turn. 
                                  % This is helpful to diagnose problems with the optimizer.
           mh_replic=0)           % Number of replications for Metropolis-Hastings algorithm                                  
           dy;                    % variables that come after the estimation command (before the semicolon) will be smoothed and plotted

% Reuse just computed mode and Hessian at the mode for initial distribution
% of Metropolis-Hastings algorithm  
estimation(datafile='simdat.mat',
           first_obs=1,           % index of first period to use
           nobs=@{nperiods},      % sample size        
           mode_compute=4,        % Choose optimizer (good choices can be 2,4,9,10)
           mode_file = 'Bayesian_RBC_leisure_mode', % This is where Dynare stores the just computed mode
           mh_replic=5000,        % Number of replications for Metropolis-Hastings algorithm                                  
           mh_nblocks = 2,        % Number of parallel chains for Metropolis-Hastings algorithm   
           mh_jscale=1.5)         % scale parameter of the jumping distribution’s covariance matrix for the Metropolis-Hastings algorithm
           dy;                    % variables that come after the estimation command (before the semicolon) will be smoothed and plotted
                                