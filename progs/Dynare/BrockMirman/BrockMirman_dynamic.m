function [residual, g1, g2, g3] = BrockMirman_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(3, 1);
T15 = params(1)*params(2)*y(6)^(-1);
T17 = T15*y(7);
T20 = y(4)^(params(1)-1);
T25 = y(1)^params(1);
lhs =y(3)^(-1);
rhs =T17*T20;
residual(1)= lhs-rhs;
lhs =y(4);
rhs =y(5)*T25-y(3);
residual(2)= lhs-rhs;
lhs =log(y(5));
rhs =params(3)*log(y(2))+params(4)*x(it_, 1);
residual(3)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(3, 8);

  %
  % Jacobian matrix
  %

  g1(1,3)=getPowerDeriv(y(3),(-1),1);
  g1(1,6)=(-(T20*y(7)*params(1)*params(2)*getPowerDeriv(y(6),(-1),1)));
  g1(1,4)=(-(T17*getPowerDeriv(y(4),params(1)-1,1)));
  g1(1,7)=(-(T15*T20));
  g1(2,3)=1;
  g1(2,1)=(-(y(5)*getPowerDeriv(y(1),params(1),1)));
  g1(2,4)=1;
  g1(2,5)=(-T25);
  g1(3,2)=(-(params(3)*1/y(2)));
  g1(3,5)=1/y(5);
  g1(3,8)=(-params(4));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],3,64);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],3,512);
end
end
end
end
