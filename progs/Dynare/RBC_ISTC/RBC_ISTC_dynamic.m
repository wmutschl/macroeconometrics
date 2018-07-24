function [residual, g1, g2, g3] = RBC_ISTC_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(9, 1);
T37 = (y(1)/y(8))^params(1);
T43 = (y(8)/y(1))^(1-params(1));
T47 = y(1)^params(1);
T48 = y(11)*T47;
T49 = y(8)^(1-params(1));
lhs =y(13)/y(5);
rhs =params(2)*y(12)/y(15)*(1-params(3)+y(14));
residual(1)= lhs-rhs;
lhs =y(5);
rhs =params(4)/(1-params(4))*(1-y(8))*y(9);
residual(2)= lhs-rhs;
lhs =y(9);
rhs =(1-params(1))*y(11)*T37;
residual(3)= lhs-rhs;
lhs =y(10);
rhs =params(1)*y(11)*T43;
residual(4)= lhs-rhs;
lhs =y(4);
rhs =T48*T49;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =(1-params(3))*y(1)+y(6);
residual(6)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(6);
residual(7)= lhs-rhs;
lhs =log(y(11));
rhs =params(5)*log(y(2))+x(it_, 1);
residual(8)= lhs-rhs;
lhs =log(y(12));
rhs =params(6)*log(y(3))+x(it_, 2);
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 17);

  %
  % Jacobian matrix
  %

T80 = getPowerDeriv(y(1)/y(8),params(1),1);
T87 = getPowerDeriv(y(8)/y(1),1-params(1),1);
  g1(1,5)=(-y(13))/(y(5)*y(5));
  g1(1,13)=1/y(5);
  g1(1,14)=(-(params(2)*y(12)/y(15)));
  g1(1,12)=(-((1-params(3)+y(14))*params(2)/y(15)));
  g1(1,15)=(-((1-params(3)+y(14))*(-(params(2)*y(12)))/(y(15)*y(15))));
  g1(2,5)=1;
  g1(2,8)=(-(y(9)*(-(params(4)/(1-params(4))))));
  g1(2,9)=(-(params(4)/(1-params(4))*(1-y(8))));
  g1(3,1)=(-((1-params(1))*y(11)*1/y(8)*T80));
  g1(3,8)=(-((1-params(1))*y(11)*T80*(-y(1))/(y(8)*y(8))));
  g1(3,9)=1;
  g1(3,11)=(-((1-params(1))*T37));
  g1(4,1)=(-(params(1)*y(11)*(-y(8))/(y(1)*y(1))*T87));
  g1(4,8)=(-(params(1)*y(11)*T87*1/y(1)));
  g1(4,10)=1;
  g1(4,11)=(-(params(1)*T43));
  g1(5,4)=1;
  g1(5,1)=(-(T49*y(11)*getPowerDeriv(y(1),params(1),1)));
  g1(5,8)=(-(T48*getPowerDeriv(y(8),1-params(1),1)));
  g1(5,11)=(-(T47*T49));
  g1(6,6)=(-1);
  g1(6,1)=(-(1-params(3)));
  g1(6,7)=1;
  g1(7,4)=1;
  g1(7,5)=(-1);
  g1(7,6)=(-1);
  g1(8,2)=(-(params(5)*1/y(2)));
  g1(8,11)=1/y(11);
  g1(8,16)=(-1);
  g1(9,3)=(-(params(6)*1/y(3)));
  g1(9,12)=1/y(12);
  g1(9,17)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,289);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,4913);
end
end
end
end
