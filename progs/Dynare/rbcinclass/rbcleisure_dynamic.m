function [residual, g1, g2, g3] = RBCleisure_dynamic(y, x, params, steady_state, it_)
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
T13 = params(4)*y(5)^(-params(8));
T18 = params(4)*params(2)*y(13)^(-params(8));
T34 = (-((-params(5))*(1-y(7))^(-params(9))));
T49 = y(2)^params(1);
T50 = y(8)*T49;
T51 = y(7)^(1-params(1));
lhs =T13;
rhs =T18*(1-params(3)+y(14));
residual(1)= lhs-rhs;
lhs =y(10);
rhs =T34/T13;
residual(2)= lhs-rhs;
lhs =y(9);
rhs =params(1)*y(4)/y(2);
residual(3)= lhs-rhs;
lhs =y(10);
rhs =y(4)*(1-params(1))/y(7);
residual(4)= lhs-rhs;
lhs =y(4);
rhs =T50*T51;
residual(5)= lhs-rhs;
lhs =y(6);
rhs =(1-params(3))*y(2)+y(11);
residual(6)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(11);
residual(7)= lhs-rhs;
lhs =log(y(8));
rhs =params(6)*log(y(3))+params(7)*x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =(y(4)-y(1))/y(1);
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 15);

  %
  % Jacobian matrix
  %

T88 = params(4)*getPowerDeriv(y(5),(-params(8)),1);
  g1(1,5)=T88;
  g1(1,13)=(-((1-params(3)+y(14))*params(4)*params(2)*getPowerDeriv(y(13),(-params(8)),1)));
  g1(1,14)=(-T18);
  g1(2,5)=(-((-(T34*T88))/(T13*T13)));
  g1(2,7)=(-((-((-params(5))*(-(getPowerDeriv(1-y(7),(-params(9)),1)))))/T13));
  g1(2,10)=1;
  g1(3,4)=(-(params(1)/y(2)));
  g1(3,2)=(-((-(params(1)*y(4)))/(y(2)*y(2))));
  g1(3,9)=1;
  g1(4,4)=(-((1-params(1))/y(7)));
  g1(4,7)=(-((-(y(4)*(1-params(1))))/(y(7)*y(7))));
  g1(4,10)=1;
  g1(5,4)=1;
  g1(5,2)=(-(T51*y(8)*getPowerDeriv(y(2),params(1),1)));
  g1(5,7)=(-(T50*getPowerDeriv(y(7),1-params(1),1)));
  g1(5,8)=(-(T49*T51));
  g1(6,2)=(-(1-params(3)));
  g1(6,6)=1;
  g1(6,11)=(-1);
  g1(7,4)=1;
  g1(7,5)=(-1);
  g1(7,11)=(-1);
  g1(8,3)=(-(params(6)*1/y(3)));
  g1(8,8)=1/y(8);
  g1(8,15)=(-params(7));
  g1(9,1)=(-(((-y(1))-(y(4)-y(1)))/(y(1)*y(1))));
  g1(9,4)=(-(1/y(1)));
  g1(9,12)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,225);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,3375);
end
end
end
end
