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
M_.fname = 'steady1_RBC_leisure';
M_.dynare_version = '4.5.3';
oo_.dynare_version = '4.5.3';
options_.dynare_version = '4.5.3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('steady1_RBC_leisure.log');
M_.exo_names = 'eps_A';
M_.exo_names_tex = '{\varepsilon_A}';
M_.exo_names_long = 'TFP shock';
M_.endo_names = 'Y';
M_.endo_names_tex = '{Y}';
M_.endo_names_long = 'output';
M_.endo_names = char(M_.endo_names, 'C');
M_.endo_names_tex = char(M_.endo_names_tex, '{C}');
M_.endo_names_long = char(M_.endo_names_long, 'consumption');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, '{K}');
M_.endo_names_long = char(M_.endo_names_long, 'capital');
M_.endo_names = char(M_.endo_names, 'L');
M_.endo_names_tex = char(M_.endo_names_tex, '{L}');
M_.endo_names_long = char(M_.endo_names_long, 'hours');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, '{Z}');
M_.endo_names_long = char(M_.endo_names_long, 'total factor productivity');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, '{R}');
M_.endo_names_long = char(M_.endo_names_long, 'interest rate');
M_.endo_names = char(M_.endo_names, 'W');
M_.endo_names_tex = char(M_.endo_names_tex, '{W}');
M_.endo_names_long = char(M_.endo_names_long, 'real wage');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, '{I}');
M_.endo_names_long = char(M_.endo_names_long, 'investment');
M_.endo_partitions = struct();
M_.param_names = 'alph';
M_.param_names_tex = '{\alpha}';
M_.param_names_long = 'capital share';
M_.param_names = char(M_.param_names, 'betta');
M_.param_names_tex = char(M_.param_names_tex, '{\beta}');
M_.param_names_long = char(M_.param_names_long, 'discount factor');
M_.param_names = char(M_.param_names, 'delt');
M_.param_names_tex = char(M_.param_names_tex, '{\delta}');
M_.param_names_long = char(M_.param_names_long, 'depreciation rate');
M_.param_names = char(M_.param_names, 'gam');
M_.param_names_tex = char(M_.param_names_tex, '{\gamma}');
M_.param_names_long = char(M_.param_names_long, 'consumption utility parameter');
M_.param_names = char(M_.param_names, 'pssi');
M_.param_names_tex = char(M_.param_names_tex, '{\psi}');
M_.param_names_long = char(M_.param_names_long, 'labor disutility parameter');
M_.param_names = char(M_.param_names, 'rhoA');
M_.param_names_tex = char(M_.param_names_tex, '{\rho_A}');
M_.param_names_long = char(M_.param_names_long, 'persistence TFP');
M_.param_names = char(M_.param_names, 'sigA');
M_.param_names_tex = char(M_.param_names_tex, '{\sigma_A}');
M_.param_names_long = char(M_.param_names_long, 'standard deviation TFP shock');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 8;
M_.param_nbr = 7;
M_.orig_endo_nbr = 8;
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
M_.hessian_eq_zero = 0;
erase_compiled_function('steady1_RBC_leisure_static');
erase_compiled_function('steady1_RBC_leisure_dynamic');
M_.orig_eq_nbr = 8;
M_.eq_nbr = 8;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 0 3 0;
 0 4 11;
 1 5 0;
 0 6 0;
 2 7 0;
 0 8 0;
 0 9 0;
 0 10 0;]';
M_.nstatic = 5;
M_.nfwrd   = 1;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 2;
M_.ndynamic   = 3;
M_.equations_tags = {
  1 , 'name' , 'Euler equation' ;
  2 , 'name' , 'labor supply' ;
  3 , 'name' , 'capital law of motion' ;
  4 , 'name' , 'market clearing/resource constraint' ;
  5 , 'name' , 'production function' ;
  6 , 'name' , 'labor demand' ;
  7 , 'name' , 'capital demand' ;
  8 , 'name' , 'exogenous TFP process' ;
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(8, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(7, 1);
M_.NNZDerivatives = [25; 24; -1];
K_o_Y = 10;   
I_o_Y = 0.25; 
M_.params( 1 ) = 0.35;
alph = M_.params( 1 );
M_.params( 4 ) = 1;
gam = M_.params( 4 );
Lbar  = 1/3;  
M_.params( 6 ) = 0.9;
rhoA = M_.params( 6 );
M_.params( 7 ) = 0.6;
sigA = M_.params( 7 );
M_.params( 3 ) = I_o_Y/K_o_Y;
delt = M_.params( 3 );
M_.params( 2 ) = 1/(1+M_.params(1)/K_o_Y-M_.params(3));
betta = M_.params( 2 );
Abar = 1;                                  
Rbar = 1/betta+delt-1;                     
K_o_L = (alph*Abar/Rbar)^(1/(1-alph));     
Kbar = K_o_L*Lbar;                         
Ybar = Kbar/K_o_Y;                         
Ibar = delt*Kbar;                          
Ybar = Abar*Kbar^alph*Lbar^(1-alph);       
Wbar = (1-alph)*Abar*(K_o_L)^alph;         
Cbar = Ybar - Ibar;                        
M_.params( 5 ) = M_.params(4)*(Cbar/Lbar)^(-1)*Wbar*(Lbar/(1-Lbar))^(-1);
pssi = M_.params( 5 );
steady;
resid;  
save('steady1_RBC_leisure_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('steady1_RBC_leisure_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
