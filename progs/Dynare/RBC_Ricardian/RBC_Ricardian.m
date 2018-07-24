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
M_.fname = 'RBC_Ricardian';
M_.dynare_version = '4.5.3';
oo_.dynare_version = '4.5.3';
options_.dynare_version = '4.5.3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('RBC_Ricardian.log');
M_.exo_names = 'epsA';
M_.exo_names_tex = 'epsA';
M_.exo_names_long = 'epsA';
M_.endo_names = 'Y';
M_.endo_names_tex = 'Y';
M_.endo_names_long = 'Y';
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, 'C');
M_.endo_names_long = char(M_.endo_names_long, 'C');
M_.endo_names = char(M_.endo_names, 'Ci');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ci');
M_.endo_names_long = char(M_.endo_names_long, 'Ci');
M_.endo_names = char(M_.endo_names, 'Cj');
M_.endo_names_tex = char(M_.endo_names_tex, 'Cj');
M_.endo_names_long = char(M_.endo_names_long, 'Cj');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'Ii');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ii');
M_.endo_names_long = char(M_.endo_names_long, 'Ii');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'Ki');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ki');
M_.endo_names_long = char(M_.endo_names_long, 'Ki');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, 'L');
M_.endo_names_long = char(M_.endo_names_long, 'L');
M_.endo_names = char(M_.endo_names, 'Li');
M_.endo_names_tex = char(M_.endo_names_tex, 'Li');
M_.endo_names_long = char(M_.endo_names_long, 'Li');
M_.endo_names = char(M_.endo_names, 'Lj');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lj');
M_.endo_names_long = char(M_.endo_names_long, 'Lj');
M_.endo_names = char(M_.endo_names, 'W');
M_.endo_names_tex = char(M_.endo_names_tex, 'W');
M_.endo_names_long = char(M_.endo_names_long, 'W');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
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
M_.param_names = char(M_.param_names, 'omeg');
M_.param_names_tex = char(M_.param_names_tex, 'omeg');
M_.param_names_long = char(M_.param_names_long, 'omeg');
M_.param_names = char(M_.param_names, 'rhoA');
M_.param_names_tex = char(M_.param_names_tex, 'rhoA');
M_.param_names_long = char(M_.param_names_long, 'rhoA');
M_.param_names = char(M_.param_names, 'sigA');
M_.param_names_tex = char(M_.param_names_tex, 'sigA');
M_.param_names_long = char(M_.param_names_long, 'sigA');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 14;
M_.param_nbr = 7;
M_.orig_endo_nbr = 14;
M_.aux_vars = [];
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
erase_compiled_function('RBC_Ricardian_static');
erase_compiled_function('RBC_Ricardian_dynamic');
M_.orig_eq_nbr = 14;
M_.eq_nbr = 14;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 0 4 0;
 0 5 0;
 0 6 18;
 0 7 0;
 0 8 0;
 0 9 0;
 1 10 0;
 2 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 19;
 3 17 0;]';
M_.nstatic = 9;
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
oo_.steady_state = zeros(14, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.NNZDerivatives = [43; -1; -1];
M_.params( 1 ) = 0.35;
alph = M_.params( 1 );
M_.params( 2 ) = 0.97;
betta = M_.params( 2 );
M_.params( 3 ) = 0.06;
delt = M_.params( 3 );
M_.params( 4 ) = 0.4;
gam = M_.params( 4 );
M_.params( 5 ) = 1;
omeg = M_.params( 5 );
M_.params( 6 ) = 0.95;
rhoA = M_.params( 6 );
M_.params( 7 ) = 0.01;
sigA = M_.params( 7 );
phhi = 0;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 1 ) = 1;
oo_.steady_state( 2 ) = 0.8;
oo_.steady_state( 3 ) = 0.6;
oo_.steady_state( 4 ) = 0.2;
oo_.steady_state( 9 ) = 0.3;
oo_.steady_state( 10 ) = 0.3;
oo_.steady_state( 11 ) = 0.3;
oo_.steady_state( 7 ) = 3.5;
oo_.steady_state( 5 ) = 0.2;
oo_.steady_state( 6 ) = 0.3;
oo_.steady_state( 12 ) = (1-M_.params(1))*oo_.steady_state(1)/oo_.steady_state(9);
oo_.steady_state( 13 ) = M_.params(1)*oo_.steady_state(1)/oo_.steady_state(7);
oo_.steady_state( 14 ) = 1;
oo_.exo_steady_state( 1 ) = 0;
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = M_.params(7)^2;
options_.order = 1;
var_list_ = char('Y','C','Ci','Cj','I','K','L','Li','Lj','W');
info = stoch_simul(var_list_);
save('RBC_Ricardian_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('RBC_Ricardian_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('RBC_Ricardian_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('RBC_Ricardian_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('RBC_Ricardian_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('RBC_Ricardian_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('RBC_Ricardian_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
