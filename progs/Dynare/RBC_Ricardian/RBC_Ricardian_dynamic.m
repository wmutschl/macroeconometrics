function [residual, g1, g2, g3] = RBC_Ricardian_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(14, 1);
T61 = (y(1)/y(12))^params(1);
T67 = (y(12)/y(1))^(1-params(1));
T71 = y(1)^params(1);
T72 = y(17)*T71;
T73 = y(12)^(1-params(1));
lhs =y(18)/y(6);
rhs =params(2)*(1-params(3)+y(19));
residual(1)= lhs-rhs;
lhs =y(7);
rhs =y(15)*y(14);
residual(2)= lhs-rhs;
lhs =y(6);
rhs =y(15)*params(4)/(1-params(4))*(1-y(13));
residual(3)= lhs-rhs;
lhs =y(7);
rhs =y(15)*params(4)/(1-params(4))*(1-y(14));
residual(4)= lhs-rhs;
lhs =y(12);
rhs =y(13)*params(5)+y(14)*(1-params(5));
residual(5)= lhs-rhs;
lhs =y(5);
rhs =y(6)*params(5)+y(7)*(1-params(5));
residual(6)= lhs-rhs;
lhs =y(10);
rhs =params(5)*y(11);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =params(5)*y(9);
residual(8)= lhs-rhs;
lhs =y(15);
rhs =(1-params(1))*y(17)*T61;
residual(9)= lhs-rhs;
lhs =y(16);
rhs =params(1)*y(17)*T67;
residual(10)= lhs-rhs;
lhs =y(4);
rhs =T72*T73;
residual(11)= lhs-rhs;
lhs =y(11);
rhs =y(9)+(1-params(3))*y(2);
residual(12)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(8);
residual(13)= lhs-rhs;
lhs =log(y(17));
rhs =params(6)*log(y(3))+x(it_, 1);
residual(14)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(14, 20);

  %
  % Jacobian matrix
  %

T97 = getPowerDeriv(y(1)/y(12),params(1),1);
T104 = getPowerDeriv(y(12)/y(1),1-params(1),1);
  g1(1,6)=(-y(18))/(y(6)*y(6));
  g1(1,18)=1/y(6);
  g1(1,19)=(-params(2));
  g1(2,7)=1;
  g1(2,14)=(-y(15));
  g1(2,15)=(-y(14));
  g1(3,6)=1;
  g1(3,13)=(-(y(15)*(-(params(4)/(1-params(4))))));
  g1(3,15)=(-(params(4)/(1-params(4))*(1-y(13))));
  g1(4,7)=1;
  g1(4,14)=(-(y(15)*(-(params(4)/(1-params(4))))));
  g1(4,15)=(-(params(4)/(1-params(4))*(1-y(14))));
  g1(5,12)=1;
  g1(5,13)=(-params(5));
  g1(5,14)=(-(1-params(5)));
  g1(6,5)=1;
  g1(6,6)=(-params(5));
  g1(6,7)=(-(1-params(5)));
  g1(7,10)=1;
  g1(7,11)=(-params(5));
  g1(8,8)=1;
  g1(8,9)=(-params(5));
  g1(9,1)=(-((1-params(1))*y(17)*1/y(12)*T97));
  g1(9,12)=(-((1-params(1))*y(17)*T97*(-y(1))/(y(12)*y(12))));
  g1(9,15)=1;
  g1(9,17)=(-((1-params(1))*T61));
  g1(10,1)=(-(params(1)*y(17)*(-y(12))/(y(1)*y(1))*T104));
  g1(10,12)=(-(params(1)*y(17)*T104*1/y(1)));
  g1(10,16)=1;
  g1(10,17)=(-(params(1)*T67));
  g1(11,4)=1;
  g1(11,1)=(-(T73*y(17)*getPowerDeriv(y(1),params(1),1)));
  g1(11,12)=(-(T72*getPowerDeriv(y(12),1-params(1),1)));
  g1(11,17)=(-(T71*T73));
  g1(12,9)=(-1);
  g1(12,2)=(-(1-params(3)));
  g1(12,11)=1;
  g1(13,4)=1;
  g1(13,5)=(-1);
  g1(13,8)=(-1);
  g1(14,3)=(-(params(6)*1/y(3)));
  g1(14,17)=1/y(17);
  g1(14,20)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,400);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,8000);
end
end
end
end
