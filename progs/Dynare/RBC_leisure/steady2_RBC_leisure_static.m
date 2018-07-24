function [residual, g1, g2, g3] = steady2_RBC_leisure_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 8, 1);

%
% Model equations
%

UC__ = params(4)*y(2)^(-params(8));
UCp__ = params(4)*y(2)^(-params(8));
UL__ = (-params(5))*(1-y(4))^(-params(9));
T47 = y(3)^params(1);
T48 = y(5)*T47;
T50 = y(4)^(1-params(1));
lhs =UC__;
rhs =params(2)*UCp__*(1-params(3)+y(6));
residual(1)= lhs-rhs;
lhs =y(7);
rhs =(-UL__)/UC__;
residual(2)= lhs-rhs;
lhs =y(3);
rhs =(1-params(3))*y(3)+y(8);
residual(3)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(8);
residual(4)= lhs-rhs;
lhs =y(1);
rhs =T48*T50;
residual(5)= lhs-rhs;
lhs =y(7);
rhs =y(1)*(1-params(1))/y(4);
residual(6)= lhs-rhs;
lhs =y(6);
rhs =y(1)*params(1)/y(3);
residual(7)= lhs-rhs;
lhs =log(y(5));
rhs =log(y(5))*params(6)+x(1);
residual(8)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(8, 8);

  %
  % Jacobian matrix
  %

T70 = params(4)*getPowerDeriv(y(2),(-params(8)),1);
  g1(1,2)=T70-(1-params(3)+y(6))*params(2)*T70;
  g1(1,6)=(-(params(2)*UCp__));
  g1(2,2)=(-((-((-UL__)*T70))/(UC__*UC__)));
  g1(2,4)=(-((-((-params(5))*(-(getPowerDeriv(1-y(4),(-params(9)),1)))))/UC__));
  g1(2,7)=1;
  g1(3,3)=1-(1-params(3));
  g1(3,8)=(-1);
  g1(4,1)=1;
  g1(4,2)=(-1);
  g1(4,8)=(-1);
  g1(5,1)=1;
  g1(5,3)=(-(T50*y(5)*getPowerDeriv(y(3),params(1),1)));
  g1(5,4)=(-(T48*getPowerDeriv(y(4),1-params(1),1)));
  g1(5,5)=(-(T47*T50));
  g1(6,1)=(-((1-params(1))/y(4)));
  g1(6,4)=(-((-(y(1)*(1-params(1))))/(y(4)*y(4))));
  g1(6,7)=1;
  g1(7,1)=(-(params(1)/y(3)));
  g1(7,3)=(-((-(y(1)*params(1)))/(y(3)*y(3))));
  g1(7,6)=1;
  g1(8,5)=1/y(5)-params(6)*1/y(5);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],8,64);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],8,512);
end
end
end
end
