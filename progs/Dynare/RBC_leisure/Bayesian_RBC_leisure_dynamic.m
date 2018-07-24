function [residual, g1, g2, g3] = Bayesian_RBC_leisure_dynamic(y, x, params, steady_state, it_)
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
UC__ = params(4)*y(5)^(-1);
UCp__ = params(4)*y(13)^(-1);
UL__ = (-params(5))*(1-y(7))^(-1);
T47 = y(2)^params(1);
T48 = y(8)*T47;
T50 = y(7)^(1-params(1));
lhs =UC__;
rhs =params(2)*UCp__*(1-params(3)+y(9));
residual(1)= lhs-rhs;
lhs =y(10);
rhs =(-UL__)/UC__;
residual(2)= lhs-rhs;
lhs =y(6);
rhs =(1-params(3))*y(2)+y(11);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(5)+y(11);
residual(4)= lhs-rhs;
lhs =y(4);
rhs =T48*T50;
residual(5)= lhs-rhs;
lhs =y(10);
rhs =y(4)*(1-params(1))/y(7);
residual(6)= lhs-rhs;
lhs =y(9);
rhs =y(4)*params(1)/y(2);
residual(7)= lhs-rhs;
lhs =log(y(8));
rhs =params(6)*log(y(3))+x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =y(4)/y(1)-1;
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 14);

  %
  % Jacobian matrix
  %

T83 = params(4)*getPowerDeriv(y(5),(-1),1);
  g1(1,5)=T83;
  g1(1,13)=(-((1-params(3)+y(9))*params(2)*params(4)*getPowerDeriv(y(13),(-1),1)));
  g1(1,9)=(-(params(2)*UCp__));
  g1(2,5)=(-((-((-UL__)*T83))/(UC__*UC__)));
  g1(2,7)=(-((-((-params(5))*(-(getPowerDeriv(1-y(7),(-1),1)))))/UC__));
  g1(2,10)=1;
  g1(3,2)=(-(1-params(3)));
  g1(3,6)=1;
  g1(3,11)=(-1);
  g1(4,4)=1;
  g1(4,5)=(-1);
  g1(4,11)=(-1);
  g1(5,4)=1;
  g1(5,2)=(-(T50*y(8)*getPowerDeriv(y(2),params(1),1)));
  g1(5,7)=(-(T48*getPowerDeriv(y(7),1-params(1),1)));
  g1(5,8)=(-(T47*T50));
  g1(6,4)=(-((1-params(1))/y(7)));
  g1(6,7)=(-((-(y(4)*(1-params(1))))/(y(7)*y(7))));
  g1(6,10)=1;
  g1(7,4)=(-(params(1)/y(2)));
  g1(7,2)=(-((-(y(4)*params(1)))/(y(2)*y(2))));
  g1(7,9)=1;
  g1(8,3)=(-(params(6)*1/y(3)));
  g1(8,8)=1/y(8);
  g1(8,14)=(-1);
  g1(9,1)=(-((-y(4))/(y(1)*y(1))));
  g1(9,4)=(-(1/y(1)));
  g1(9,12)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,2744);
end
end
end
end
