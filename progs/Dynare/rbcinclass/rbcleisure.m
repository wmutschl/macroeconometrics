%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'RBCleisure';
M_.dynare_version = '4.5.3';
oo_.dynare_version = '4.5.3';
options_.dynare_version = '4.5.3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('RBCleisure.log');
M_.exo_names = 'eps_A';
M_.exo_names_tex = 'eps\_A';
M_.exo_names_long = 'eps_A';
M_.endo_names = 'Y';
M_.endo_names_tex = 'Y';
M_.endo_names_long = 'Y';
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, 'C');
M_.endo_names_long = char(M_.endo_names_long, 'C');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'W');
M_.endo_names_tex = char(M_.endo_names_tex, 'W');
M_.endo_names_long = char(M_.endo_names_long, 'W');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'growth');
M_.endo_names_tex = char(M_.endo_names_tex, 'growth');
M_.endo_names_long = char(M_.endo_names_long, 'growth');
M_.endo_partitions = struct();
M_.param_names = 'alph';
M_.param_names_tex = 'alph';
M_.param_names_long = 'alph';
M_.param_names = char(M_.param_names, 'betta');
M_.param_names_tex = char(M_.param_names_tex, 'betta');
M_.param_names_long = char(M_.param_names_long, 'betta');
M_.param_names = char(M_.param_names, 'delt');
M_.param_names_tex = char(M_.param_names_tex, 'delt');
M_.param_names_long = char(M_.param_names_long, 'delt');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, 'gam');
M_.param_names_long = char(M_.param_names_long, 'gam');
M_.param_names = char(M_.param_names, 'pssi');
M_.param_names_tex = char(M_.param_names_tex, 'pssi');
M_.param_names_long = char(M_.param_names_long, 'pssi');
M_.param_names = char(M_.param_names, 'rhoA');
M_.param_names_tex = char(M_.param_names_tex, 'rhoA');
M_.param_names_long = char(M_.param_names_long, 'rhoA');
M_.param_names = char(M_.param_names, 'sigA');
M_.param_names_tex = char(M_.param_names_tex, 'sigA');
M_.param_names_long = char(M_.param_names_long, 'sigA');
M_.param_names = char(M_.param_names, 'etaC');
M_.param_names_tex = char(M_.param_names_tex, 'etaC');
M_.param_names_long = char(M_.param_names_long, 'etaC');
M_.param_names = char(M_.param_names, 'etaL');
M_.param_names_tex = char(M_.param_names_tex, 'etaL');
M_.param_names_long = char(M_.param_names_long, 'etaL');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 9;
M_.param_nbr = 9;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'growth'};
options_.varobs_id = [ 9  ];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('RBCleisure_static');
erase_compiled_function('RBCleisure_dynamic');
M_.orig_eq_nbr = 9;
M_.eq_nbr = 9;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 1 4 0;
 0 5 13;
 2 6 0;
 0 7 0;
 3 8 0;
 0 9 14;
 0 10 0;
 0 11 0;
 0 12 0;]';
M_.nstatic = 4;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(9, 1);
M_.NNZDerivatives = [28; -1; -1];
K_o_Y = 10; 
I_o_Y = 0.25; 
M_.params( 1 ) = 0.35;
alph = M_.params( 1 );
M_.params( 4 ) = 1;
gam = M_.params( 4 );
Lbar = 1/3;
M_.params( 6 ) = 0.9;
rhoA = M_.params( 6 );
M_.params( 7 ) = 0.6;
sigA = M_.params( 7 );
M_.params( 8 ) = 1;
etaC = M_.params( 8 );
M_.params( 9 ) = 1;
etaL = M_.params( 9 );
M_.params( 3 ) = I_o_Y/K_o_Y;
delt = M_.params( 3 );
M_.params( 2 ) = 1/(1+M_.params(1)/K_o_Y-M_.params(3));
betta = M_.params( 2 );
Abar = 1;
Rbar = 1/betta+delt-1;
K_o_L = ((alph*Abar)/Rbar)^(1/(1-alph));
Kbar = K_o_L*Lbar;
Ybar = Kbar/K_o_Y;
Ibar = delt*Kbar;
Wbar = (1-alph)*Abar*K_o_L^alph;
Cbar = Ybar - Ibar;
M_.params( 5 ) = M_.params(4)*(Cbar/Lbar)^(-M_.params(8))*Wbar*(Lbar/(1-Lbar))^(-M_.params(9));
pssi = M_.params( 5 );
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1;
resid;
options_.solve_algo = 4;
steady;
resid;
options_.irf = 0;
options_.order = 1;
options_.periods = 10000;
var_list_ = char();
info = stoch_simul(var_list_);
save('simdat','growth');
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.param_vals = [estim_params_.param_vals; 1, NaN, (-Inf), Inf, 3, .3, .1, 0.1, 0.9, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, NaN, (-Inf), Inf, 1, 0.8, 0.1, 0.01, 0.99, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, NaN, (-Inf), Inf, 4, 0.5, 4, 0.0001, 5, NaN ];
options_.mh_replic = 0;
options_.mode_check.status = 1;
options_.mode_compute = 4;
options_.datafile = 'simdat.mat';
options_.first_obs = 203;
options_.nobs = 200;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
options_.mh_jscale = 1.5;
options_.mh_nblck = 2;
options_.mh_replic = 5000;
options_.mode_check.status = 1;
options_.mode_compute = 0;
options_.datafile = 'simdat.mat';
options_.mode_file = 'RBCleisure_mode';
options_.first_obs = 203;
options_.nobs = 200;
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
save('RBCleisure_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('RBCleisure_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('RBCleisure_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('RBCleisure_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('RBCleisure_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('RBCleisure_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('RBCleisure_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
