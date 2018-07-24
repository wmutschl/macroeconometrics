function [residual, g1, g2, g3] = simdat_RBC_leisure_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(8, 1);
UC__ = params(4)*y(4)^(-1);
UCp__ = params(4)*y(11)^(-1);
UL__ = (-params(5))*(1-y(6))^(-1);
T47 = y(1)^params(1);
T48 = y(7)*T47;
T50 = y(6)^(1-params(1));
lhs =UC__;
rhs =params(2)*UCp__*(1-params(3)+y(8));
residual(1)= lhs-rhs;
lhs =y(9);
rhs =(-UL__)/UC__;
residual(2)= lhs-rhs;
lhs =y(5);
rhs =(1-params(3))*y(1)+y(10);
residual(3)= lhs-rhs;
lhs =y(3);
rhs =y(4)+y(10);
residual(4)= lhs-rhs;
lhs =y(3);
rhs =T48*T50;
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(3)*(1-params(1))/y(6);
residual(6)= lhs-rhs;
lhs =y(8);
rhs =y(3)*params(1)/y(1);
residual(7)= lhs-rhs;
lhs =log(y(7));
rhs =params(6)*log(y(2))+x(it_, 1);
residual(8)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(8, 12);

  %
  % Jacobian matrix
  %

T72 = params(4)*getPowerDeriv(y(4),(-1),1);
  g1(1,4)=T72;
  g1(1,11)=(-((1-params(3)+y(8))*params(2)*params(4)*getPowerDeriv(y(11),(-1),1)));
  g1(1,8)=(-(params(2)*UCp__));
  g1(2,4)=(-((-((-UL__)*T72))/(UC__*UC__)));
  g1(2,6)=(-((-((-params(5))*(-(getPowerDeriv(1-y(6),(-1),1)))))/UC__));
  g1(2,9)=1;
  g1(3,1)=(-(1-params(3)));
  g1(3,5)=1;
  g1(3,10)=(-1);
  g1(4,3)=1;
  g1(4,4)=(-1);
  g1(4,10)=(-1);
  g1(5,3)=1;
  g1(5,1)=(-(T50*y(7)*getPowerDeriv(y(1),params(1),1)));
  g1(5,6)=(-(T48*getPowerDeriv(y(6),1-params(1),1)));
  g1(5,7)=(-(T47*T50));
  g1(6,3)=(-((1-params(1))/y(6)));
  g1(6,6)=(-((-(y(3)*(1-params(1))))/(y(6)*y(6))));
  g1(6,9)=1;
  g1(7,3)=(-(params(1)/y(1)));
  g1(7,1)=(-((-(y(3)*params(1)))/(y(1)*y(1))));
  g1(7,8)=1;
  g1(8,2)=(-(params(6)*1/y(2)));
  g1(8,7)=1/y(7);
  g1(8,12)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,144);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,1728);
end
end
end
end
