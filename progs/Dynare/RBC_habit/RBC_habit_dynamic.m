function [residual, g1, g2, g3] = RBC_habit_dynamic(y, x, params, steady_state, it_)
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
cback__ = 1/(y(5)-params(6)*y(1));
c__ = 1/(y(13)-y(5)*params(6));
cp__ = y(15);
T57 = (y(2)/y(8))^params(1);
T63 = (y(8)/y(2))^(1-params(1));
T67 = y(2)^params(1);
T68 = y(11)*T67;
T69 = y(8)^(1-params(1));
lhs =(cback__-params(6)*params(2)*c__)/(c__-params(6)*params(2)*cp__);
rhs =params(2)*(1+y(14)-params(3));
residual(1)= lhs-rhs;
lhs =(cback__*params(4)-c__*params(6)*params(2)*params(4))*y(9);
rhs =(1-params(4))/(1-y(8));
residual(2)= lhs-rhs;
lhs =y(9);
rhs =(1-params(1))*y(11)*T57;
residual(3)= lhs-rhs;
lhs =y(10);
rhs =params(1)*y(11)*T63;
residual(4)= lhs-rhs;
lhs =y(4);
rhs =T68*T69;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =y(2)*(1-params(3))+y(6);
residual(6)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(6);
residual(7)= lhs-rhs;
lhs =log(y(11));
rhs =params(5)*log(y(3))+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =1/(y(13)-y(5)*params(6));
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 16);

  %
  % Jacobian matrix
  %

T99 = params(6)/((y(13)-y(5)*params(6))*(y(13)-y(5)*params(6)));
T112 = (-1)/((y(13)-y(5)*params(6))*(y(13)-y(5)*params(6)));
T124 = getPowerDeriv(y(2)/y(8),params(1),1);
T131 = getPowerDeriv(y(8)/y(2),1-params(1),1);
  g1(1,1)=params(6)/((y(5)-params(6)*y(1))*(y(5)-params(6)*y(1)))/(c__-params(6)*params(2)*cp__);
  g1(1,5)=((c__-params(6)*params(2)*cp__)*((-1)/((y(5)-params(6)*y(1))*(y(5)-params(6)*y(1)))-params(6)*params(2)*T99)-(cback__-params(6)*params(2)*c__)*T99)/((c__-params(6)*params(2)*cp__)*(c__-params(6)*params(2)*cp__));
  g1(1,13)=((c__-params(6)*params(2)*cp__)*(-(params(6)*params(2)*T112))-(cback__-params(6)*params(2)*c__)*T112)/((c__-params(6)*params(2)*cp__)*(c__-params(6)*params(2)*cp__));
  g1(1,14)=(-params(2));
  g1(1,15)=(-((cback__-params(6)*params(2)*c__)*(-(params(6)*params(2)))))/((c__-params(6)*params(2)*cp__)*(c__-params(6)*params(2)*cp__));
  g1(2,1)=y(9)*params(4)*params(6)/((y(5)-params(6)*y(1))*(y(5)-params(6)*y(1)));
  g1(2,5)=y(9)*(params(4)*(-1)/((y(5)-params(6)*y(1))*(y(5)-params(6)*y(1)))-params(6)*params(2)*params(4)*T99);
  g1(2,13)=y(9)*(-(params(6)*params(2)*params(4)*T112));
  g1(2,8)=(-((1-params(4))/((1-y(8))*(1-y(8)))));
  g1(2,9)=cback__*params(4)-c__*params(6)*params(2)*params(4);
  g1(3,2)=(-((1-params(1))*y(11)*1/y(8)*T124));
  g1(3,8)=(-((1-params(1))*y(11)*T124*(-y(2))/(y(8)*y(8))));
  g1(3,9)=1;
  g1(3,11)=(-((1-params(1))*T57));
  g1(4,2)=(-(params(1)*y(11)*(-y(8))/(y(2)*y(2))*T131));
  g1(4,8)=(-(params(1)*y(11)*T131*1/y(2)));
  g1(4,10)=1;
  g1(4,11)=(-(params(1)*T63));
  g1(5,4)=1;
  g1(5,2)=(-(T69*y(11)*getPowerDeriv(y(2),params(1),1)));
  g1(5,8)=(-(T68*getPowerDeriv(y(8),1-params(1),1)));
  g1(5,11)=(-(T67*T69));
  g1(6,6)=(-1);
  g1(6,2)=(-(1-params(3)));
  g1(6,7)=1;
  g1(7,4)=1;
  g1(7,5)=(-1);
  g1(7,6)=(-1);
  g1(8,3)=(-(params(5)*1/y(3)));
  g1(8,11)=1/y(11);
  g1(8,16)=(-1);
  g1(9,5)=(-T99);
  g1(9,13)=(-T112);
  g1(9,12)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,256);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,4096);
end
end
end
end
