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
M_.fname = 'RBC_homeproduction';
M_.dynare_version = '4.5.3';
oo_.dynare_version = '4.5.3';
options_.dynare_version = '4.5.3';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('RBC_homeproduction.log');
M_.exo_names = 'epsA';
M_.exo_names_tex = 'epsA';
M_.exo_names_long = 'epsA';
M_.exo_names = char(M_.exo_names, 'epsB');
M_.exo_names_tex = char(M_.exo_names_tex, 'epsB');
M_.exo_names_long = char(M_.exo_names_long, 'epsB');
M_.endo_names = 'Y';
M_.endo_names_tex = 'Y';
M_.endo_names_long = 'Y';
M_.endo_names = char(M_.endo_names, 'Cm');
M_.endo_names_tex = char(M_.endo_names_tex, 'Cm');
M_.endo_names_long = char(M_.endo_names_long, 'Cm');
M_.endo_names = char(M_.endo_names, 'Ch');
M_.endo_names_tex = char(M_.endo_names_tex, 'Ch');
M_.endo_names_long = char(M_.endo_names_long, 'Ch');
M_.endo_names = char(M_.endo_names, 'I');
M_.endo_names_tex = char(M_.endo_names_tex, 'I');
M_.endo_names_long = char(M_.endo_names_long, 'I');
M_.endo_names = char(M_.endo_names, 'K');
M_.endo_names_tex = char(M_.endo_names_tex, 'K');
M_.endo_names_long = char(M_.endo_names_long, 'K');
M_.endo_names = char(M_.endo_names, 'Lm');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lm');
M_.endo_names_long = char(M_.endo_names_long, 'Lm');
M_.endo_names = char(M_.endo_names, 'Lh');
M_.endo_names_tex = char(M_.endo_names_tex, 'Lh');
M_.endo_names_long = char(M_.endo_names_long, 'Lh');
M_.endo_names = char(M_.endo_names, 'W');
M_.endo_names_tex = char(M_.endo_names_tex, 'W');
M_.endo_names_long = char(M_.endo_names_long, 'W');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'A');
M_.endo_names_tex = char(M_.endo_names_tex, 'A');
M_.endo_names_long = char(M_.endo_names_long, 'A');
M_.endo_names = char(M_.endo_names, 'B');
M_.endo_names_tex = char(M_.endo_names_tex, 'B');
M_.endo_names_long = char(M_.endo_names_long, 'B');
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
M_.param_names = char(M_.param_names, 'etta');
M_.param_names_tex = char(M_.param_names_tex, 'etta');
M_.param_names_long = char(M_.param_names_long, 'etta');
M_.param_names = char(M_.param_names, 'thet');
M_.param_names_tex = char(M_.param_names_tex, 'thet');
M_.param_names_long = char(M_.param_names_long, 'thet');
M_.param_names = char(M_.param_names, 'rhoA');
M_.param_names_tex = char(M_.param_names_tex, 'rhoA');
M_.param_names_long = char(M_.param_names_long, 'rhoA');
M_.param_names = char(M_.param_names, 'rhoB');
M_.param_names_tex = char(M_.param_names_tex, 'rhoB');
M_.param_names_long = char(M_.param_names_long, 'rhoB');
M_.param_names = char(M_.param_names, 'sigA');
M_.param_names_tex = char(M_.param_names_tex, 'sigA');
M_.param_names_long = char(M_.param_names_long, 'sigA');
M_.param_names = char(M_.param_names, 'sigB');
M_.param_names_tex = char(M_.param_names_tex, 'sigB');
M_.param_names_long = char(M_.param_names_long, 'sigB');
M_.param_names = char(M_.param_names, 'Abar');
M_.param_names_tex = char(M_.param_names_tex, 'Abar');
M_.param_names_long = char(M_.param_names_long, 'Abar');
M_.param_names = char(M_.param_names, 'Bbar');
M_.param_names_tex = char(M_.param_names_tex, 'Bbar');
M_.param_names_long = char(M_.param_names_long, 'Bbar');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 11;
M_.param_nbr = 13;
M_.orig_endo_nbr = 11;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=1;
options_.bytecode=1;
options_.use_dll=0;
M_.hessian_eq_zero = 0;
if exist('bytecode') ~= 3
  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')
end
erase_compiled_function('RBC_homeproduction_static');
erase_compiled_function('RBC_homeproduction_dynamic');
M_.orig_eq_nbr = 11;
M_.eq_nbr = 11;
M_.ramsey_eq_nbr = 0;
M_.lead_lag_incidence = [
 0 4 0;
 0 5 15;
 0 6 16;
 0 7 0;
 1 8 0;
 0 9 0;
 0 10 0;
 0 11 0;
 0 12 17;
 2 13 0;
 3 14 0;]';
M_.nstatic = 5;
M_.nfwrd   = 3;
M_.npred   = 3;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 6;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
block_structure.block(1).Simulation_Type = 3;
block_structure.block(1).maximum_lag = 1;
block_structure.block(1).maximum_lead = 0;
block_structure.block(1).maximum_endo_lag = 1;
block_structure.block(1).maximum_endo_lead = 0;
block_structure.block(1).maximum_exo_lag = 0;
block_structure.block(1).maximum_exo_lead = 0;
block_structure.block(1).maximum_exo_det_lag = 0;
block_structure.block(1).maximum_exo_det_lead = 0;
block_structure.block(1).endo_nbr = 1;
block_structure.block(1).mfs = 1;
block_structure.block(1).equation = [ 10];
block_structure.block(1).variable = [ 10];
block_structure.block(1).exo_nbr = 1;
block_structure.block(1).exogenous = [ 1];
block_structure.block(1).exogenous_det = [];
block_structure.block(1).exo_det_nbr = 0;
block_structure.block(1).other_endogenous = [];
block_structure.block(1).other_endogenous_block = [];
block_structure.block(1).tm1 = zeros(0, 3);
block_structure.block(1).other_endo_nbr = 0;
block_structure.block(1).lead_lag_incidence = [];
block_structure.block(1).lead_lag_incidence = [ block_structure.block(1).lead_lag_incidence; 1]; %lag = -1
block_structure.block(1).lead_lag_incidence = [ block_structure.block(1).lead_lag_incidence; 2]; %lag = 0
block_structure.block(1).lead_lag_incidence = [ block_structure.block(1).lead_lag_incidence;  0]; %lag = 1
block_structure.block(1).sorted_col_dr_ghx = [1 ];
block_structure.block(1).lead_lag_incidence_other = [];
block_structure.block(1).lead_lag_incidence_other = [ block_structure.block(1).lead_lag_incidence_other; ]; %lag = -1
block_structure.block(1).lead_lag_incidence_other = [ block_structure.block(1).lead_lag_incidence_other; ]; %lag = 0
block_structure.block(1).lead_lag_incidence_other = [ block_structure.block(1).lead_lag_incidence_other; ]; %lag = 1
block_structure.block(1).n_static = 0;
block_structure.block(1).n_forward = 0;
block_structure.block(1).n_backward = 1;
block_structure.block(1).n_mixed = 0;
block_structure.block(2).Simulation_Type = 3;
block_structure.block(2).maximum_lag = 1;
block_structure.block(2).maximum_lead = 0;
block_structure.block(2).maximum_endo_lag = 1;
block_structure.block(2).maximum_endo_lead = 0;
block_structure.block(2).maximum_exo_lag = 0;
block_structure.block(2).maximum_exo_lead = 0;
block_structure.block(2).maximum_exo_det_lag = 0;
block_structure.block(2).maximum_exo_det_lead = 0;
block_structure.block(2).endo_nbr = 1;
block_structure.block(2).mfs = 1;
block_structure.block(2).equation = [ 11];
block_structure.block(2).variable = [ 11];
block_structure.block(2).exo_nbr = 1;
block_structure.block(2).exogenous = [ 2];
block_structure.block(2).exogenous_det = [];
block_structure.block(2).exo_det_nbr = 0;
block_structure.block(2).other_endogenous = [];
block_structure.block(2).other_endogenous_block = [];
block_structure.block(2).tm1 = zeros(0, 3);
block_structure.block(2).other_endo_nbr = 0;
block_structure.block(2).lead_lag_incidence = [];
block_structure.block(2).lead_lag_incidence = [ block_structure.block(2).lead_lag_incidence; 1]; %lag = -1
block_structure.block(2).lead_lag_incidence = [ block_structure.block(2).lead_lag_incidence; 2]; %lag = 0
block_structure.block(2).lead_lag_incidence = [ block_structure.block(2).lead_lag_incidence;  0]; %lag = 1
block_structure.block(2).sorted_col_dr_ghx = [2 ];
block_structure.block(2).lead_lag_incidence_other = [];
block_structure.block(2).lead_lag_incidence_other = [ block_structure.block(2).lead_lag_incidence_other; ]; %lag = -1
block_structure.block(2).lead_lag_incidence_other = [ block_structure.block(2).lead_lag_incidence_other; ]; %lag = 0
block_structure.block(2).lead_lag_incidence_other = [ block_structure.block(2).lead_lag_incidence_other; ]; %lag = 1
block_structure.block(2).n_static = 0;
block_structure.block(2).n_forward = 0;
block_structure.block(2).n_backward = 1;
block_structure.block(2).n_mixed = 0;
block_structure.block(3).Simulation_Type = 6;
block_structure.block(3).maximum_lag = 0;
block_structure.block(3).maximum_lead = 0;
block_structure.block(3).maximum_endo_lag = 0;
block_structure.block(3).maximum_endo_lead = 0;
block_structure.block(3).maximum_exo_lag = 0;
block_structure.block(3).maximum_exo_lead = 0;
block_structure.block(3).maximum_exo_det_lag = 0;
block_structure.block(3).maximum_exo_det_lead = 0;
block_structure.block(3).endo_nbr = 3;
block_structure.block(3).mfs = 3;
block_structure.block(3).equation = [ 3 2 6];
block_structure.block(3).variable = [ 7 6 3];
block_structure.block(3).exo_nbr = 0;
block_structure.block(3).exogenous = [];
block_structure.block(3).exogenous_det = [];
block_structure.block(3).exo_det_nbr = 0;
block_structure.block(3).other_endogenous = [ 11];
block_structure.block(3).other_endogenous_block = [ 2];
block_structure.block(3).tm1 = zeros(1, 3);
block_structure.block(3).tm1(1, 2) = 1;
block_structure.block(3).other_endo_nbr = 1;
block_structure.block(3).lead_lag_incidence = [];
block_structure.block(3).lead_lag_incidence = [ block_structure.block(3).lead_lag_incidence;  0 0 0]; %lag = -1
block_structure.block(3).lead_lag_incidence = [ block_structure.block(3).lead_lag_incidence; 1 2 3]; %lag = 0
block_structure.block(3).lead_lag_incidence = [ block_structure.block(3).lead_lag_incidence;  0 0 0]; %lag = 1
block_structure.block(3).sorted_col_dr_ghx = [];
block_structure.block(3).lead_lag_incidence_other = [];
block_structure.block(3).lead_lag_incidence_other = [ block_structure.block(3).lead_lag_incidence_other;  0]; %lag = -1
block_structure.block(3).lead_lag_incidence_other = [ block_structure.block(3).lead_lag_incidence_other;  1]; %lag = 0
block_structure.block(3).lead_lag_incidence_other = [ block_structure.block(3).lead_lag_incidence_other;  0]; %lag = 1
block_structure.block(3).n_static = 3;
block_structure.block(3).n_forward = 0;
block_structure.block(3).n_backward = 0;
block_structure.block(3).n_mixed = 0;
block_structure.block(4).Simulation_Type = 8;
block_structure.block(4).maximum_lag = 1;
block_structure.block(4).maximum_lead = 1;
block_structure.block(4).maximum_endo_lag = 1;
block_structure.block(4).maximum_endo_lead = 1;
block_structure.block(4).maximum_exo_lag = 0;
block_structure.block(4).maximum_exo_lead = 0;
block_structure.block(4).maximum_exo_det_lag = 0;
block_structure.block(4).maximum_exo_det_lead = 0;
block_structure.block(4).endo_nbr = 5;
block_structure.block(4).mfs = 5;
block_structure.block(4).equation = [ 7 9 8 5 1];
block_structure.block(4).variable = [ 1 4 5 9 2];
block_structure.block(4).exo_nbr = 0;
block_structure.block(4).exogenous = [];
block_structure.block(4).exogenous_det = [];
block_structure.block(4).exo_det_nbr = 0;
block_structure.block(4).other_endogenous = [ 3 6 10];
block_structure.block(4).other_endogenous_block = [ 3 3 1];
block_structure.block(4).tm1 = zeros(3, 3);
block_structure.block(4).tm1(3, 1) = 1;
block_structure.block(4).other_endo_nbr = 3;
block_structure.block(4).lead_lag_incidence = [];
block_structure.block(4).lead_lag_incidence = [ block_structure.block(4).lead_lag_incidence;  0 0 1 0 0]; %lag = -1
block_structure.block(4).lead_lag_incidence = [ block_structure.block(4).lead_lag_incidence; 2 3 4 5 6]; %lag = 0
block_structure.block(4).lead_lag_incidence = [ block_structure.block(4).lead_lag_incidence;  0 0 0 7 8]; %lag = 1
block_structure.block(4).sorted_col_dr_ghx = [3 ];
block_structure.block(4).lead_lag_incidence_other = [];
block_structure.block(4).lead_lag_incidence_other = [ block_structure.block(4).lead_lag_incidence_other;  0 0 0]; %lag = -1
block_structure.block(4).lead_lag_incidence_other = [ block_structure.block(4).lead_lag_incidence_other;  1 2 3]; %lag = 0
block_structure.block(4).lead_lag_incidence_other = [ block_structure.block(4).lead_lag_incidence_other;  4 0 0]; %lag = 1
block_structure.block(4).n_static = 2;
block_structure.block(4).n_forward = 2;
block_structure.block(4).n_backward = 1;
block_structure.block(4).n_mixed = 0;
block_structure.block(5).Simulation_Type = 1;
block_structure.block(5).maximum_lag = 1;
block_structure.block(5).maximum_lead = 0;
block_structure.block(5).maximum_endo_lag = 0;
block_structure.block(5).maximum_endo_lead = 0;
block_structure.block(5).maximum_exo_lag = 0;
block_structure.block(5).maximum_exo_lead = 0;
block_structure.block(5).maximum_exo_det_lag = 0;
block_structure.block(5).maximum_exo_det_lead = 0;
block_structure.block(5).endo_nbr = 1;
block_structure.block(5).mfs = 1;
block_structure.block(5).equation = [ 4];
block_structure.block(5).variable = [ 8];
block_structure.block(5).exo_nbr = 0;
block_structure.block(5).exogenous = [];
block_structure.block(5).exogenous_det = [];
block_structure.block(5).exo_det_nbr = 0;
block_structure.block(5).other_endogenous = [ 5 6 10];
block_structure.block(5).other_endogenous_block = [ 4 3 1];
block_structure.block(5).tm1 = zeros(3, 3);
block_structure.block(5).tm1(1, 3) = 1;
block_structure.block(5).tm1(3, 1) = 1;
block_structure.block(5).other_endo_nbr = 3;
block_structure.block(5).lead_lag_incidence = [];
block_structure.block(5).lead_lag_incidence = [ block_structure.block(5).lead_lag_incidence;  0]; %lag = -1
block_structure.block(5).lead_lag_incidence = [ block_structure.block(5).lead_lag_incidence; 1]; %lag = 0
block_structure.block(5).lead_lag_incidence = [ block_structure.block(5).lead_lag_incidence;  0]; %lag = 1
block_structure.block(5).sorted_col_dr_ghx = [];
block_structure.block(5).lead_lag_incidence_other = [];
block_structure.block(5).lead_lag_incidence_other = [ block_structure.block(5).lead_lag_incidence_other;  1 0 0]; %lag = -1
block_structure.block(5).lead_lag_incidence_other = [ block_structure.block(5).lead_lag_incidence_other;  0 2 3]; %lag = 0
block_structure.block(5).lead_lag_incidence_other = [ block_structure.block(5).lead_lag_incidence_other;  0 0 0]; %lag = 1
block_structure.block(5).n_static = 1;
block_structure.block(5).n_forward = 0;
block_structure.block(5).n_backward = 0;
block_structure.block(5).n_mixed = 0;
M_.block_structure.block = block_structure.block;
M_.block_structure.variable_reordered = [ 10 11 7 6 3 1 4 5 9 2 8];
M_.block_structure.equation_reordered = [ 10 11 3 2 6 7 9 8 5 1 4];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [4 5;
5 5;
7 5;
8 5;
10 10;
11 11;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [1 2;
1 3;
2 6;
2 7;
3 3;
3 6;
3 7;
4 6;
4 8;
4 10;
5 6;
5 9;
5 10;
6 3;
6 7;
6 11;
7 1;
7 6;
7 10;
8 4;
8 5;
9 1;
9 2;
9 4;
10 10;
11 11;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [1 2;
1 3;
1 9;
];
M_.state_var = [10 11 5 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(11, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(13, 1);
M_.NNZDerivatives = [41; 67; -1];
block_structure_stat.block(1).Simulation_Type = 3;
block_structure_stat.block(1).endo_nbr = 1;
block_structure_stat.block(1).mfs = 1;
block_structure_stat.block(1).equation = [ 10];
block_structure_stat.block(1).variable = [ 10];
block_structure_stat.block(2).Simulation_Type = 3;
block_structure_stat.block(2).endo_nbr = 1;
block_structure_stat.block(2).mfs = 1;
block_structure_stat.block(2).equation = [ 11];
block_structure_stat.block(2).variable = [ 11];
block_structure_stat.block(3).Simulation_Type = 6;
block_structure_stat.block(3).endo_nbr = 3;
block_structure_stat.block(3).mfs = 1;
block_structure_stat.block(3).equation = [ 3 2 6];
block_structure_stat.block(3).variable = [ 6 7 3];
block_structure_stat.block(4).Simulation_Type = 6;
block_structure_stat.block(4).endo_nbr = 5;
block_structure_stat.block(4).mfs = 1;
block_structure_stat.block(4).equation = [ 5 7 8 9 1];
block_structure_stat.block(4).variable = [ 5 1 4 2 9];
block_structure_stat.block(5).Simulation_Type = 1;
block_structure_stat.block(5).endo_nbr = 1;
block_structure_stat.block(5).mfs = 1;
block_structure_stat.block(5).equation = [ 4];
block_structure_stat.block(5).variable = [ 8];
M_.block_structure_stat.block = block_structure_stat.block;
M_.block_structure_stat.variable_reordered = [ 10 11 6 7 3 5 1 4 2 9 8];
M_.block_structure_stat.equation_reordered = [ 10 11 3 2 6 5 7 8 9 1 4];
M_.block_structure_stat.incidence.sparse_IM = [1 2;
1 3;
1 9;
2 6;
2 7;
3 3;
3 6;
3 7;
4 5;
4 6;
4 8;
4 10;
5 5;
5 6;
5 9;
5 10;
6 3;
6 7;
6 11;
7 1;
7 5;
7 6;
7 10;
8 4;
8 5;
9 1;
9 2;
9 4;
10 10;
11 11;
];
M_.params( 1 ) = 0.35;
alph = M_.params( 1 );
M_.params( 2 ) = 0.97;
betta = M_.params( 2 );
M_.params( 3 ) = 0.06;
delt = M_.params( 3 );
M_.params( 4 ) = 0.4;
gam = M_.params( 4 );
M_.params( 5 ) = 0;
omeg = M_.params( 5 );
M_.params( 6 ) = 0.8;
etta = M_.params( 6 );
M_.params( 7 ) = 1;
thet = M_.params( 7 );
M_.params( 8 ) = 0.95;
rhoA = M_.params( 8 );
M_.params( 10 ) = 0.01;
sigA = M_.params( 10 );
M_.params( 9 ) = 0.95;
rhoB = M_.params( 9 );
M_.params( 11 ) = 0.01;
sigB = M_.params( 11 );
M_.params( 12 ) = 1;
Abar = M_.params( 12 );
M_.params( 13 ) = 1;
Bbar = M_.params( 13 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.exo_steady_state( 1 ) = 0;
oo_.exo_steady_state( 2 ) = 0;
oo_.steady_state( 10 ) = 1;
oo_.steady_state( 11 ) = 1;
oo_.steady_state( 1 ) = 0.459779;
oo_.steady_state( 2 ) = 0.353592;
oo_.steady_state( 3 ) = 0.199197;
oo_.steady_state( 6 ) = 0.22251;
oo_.steady_state( 7 ) = 0.133077;
oo_.steady_state( 5 ) = 1.76978;
oo_.steady_state( 4 ) = 0.106187;
oo_.steady_state( 8 ) = (1-M_.params(1))*oo_.steady_state(1)/oo_.steady_state(6);
oo_.steady_state( 9 ) = M_.params(1)*oo_.steady_state(1)/oo_.steady_state(5);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
steady;
save('RBC_homeproduction_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('RBC_homeproduction_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('RBC_homeproduction_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('RBC_homeproduction_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('RBC_homeproduction_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('RBC_homeproduction_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('RBC_homeproduction_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
