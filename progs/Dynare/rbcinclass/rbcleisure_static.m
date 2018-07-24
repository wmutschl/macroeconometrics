function [residual, g1, g2, g3] = RBCleisure_static(y, x, params)
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

residual = zeros( 9, 1);

%
% Model equations
%

T12 = y(2)^(-params(8));
T32 = (-((-params(5))*(1-y(4))^(-params(9))));
T46 = y(3)^params(1);
T47 = y(5)*T46;
T48 = y(4)^(1-params(1));
lhs =params(4)*T12;
rhs =T12*params(4)*params(2)*(1-params(3)+y(6));
residual(1)= lhs-rhs;
lhs =y(7);
rhs =T32/(params(4)*T12);
residual(2)= lhs-rhs;
lhs =y(6);
rhs =params(1)*y(1)/y(3);
residual(3)= lhs-rhs;
lhs =y(7);
rhs =y(1)*(1-params(1))/y(4);
residual(4)= lhs-rhs;
lhs =y(1);
rhs =T47*T48;
residual(5)= lhs-rhs;
lhs =y(3);
rhs =(1-params(3))*y(3)+y(8);
residual(6)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(8);
residual(7)= lhs-rhs;
lhs =log(y(5));
rhs =log(y(5))*params(6)+params(7)*x(1);
residual(8)= lhs-rhs;
residual(9) = y(9);
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

  %
  % Jacobian matrix
  %

T71 = getPowerDeriv(y(2),(-params(8)),1);
  g1(1,2)=params(4)*T71-(1-params(3)+y(6))*params(4)*params(2)*T71;
  g1(1,6)=(-(T12*params(4)*params(2)));
  g1(2,2)=(-((-(T32*params(4)*T71))/(params(4)*T12*params(4)*T12)));
  g1(2,4)=(-((-((-params(5))*(-(getPowerDeriv(1-y(4),(-params(9)),1)))))/(params(4)*T12)));
  g1(2,7)=1;
  g1(3,1)=(-(params(1)/y(3)));
  g1(3,3)=(-((-(params(1)*y(1)))/(y(3)*y(3))));
  g1(3,6)=1;
  g1(4,1)=(-((1-params(1))/y(4)));
  g1(4,4)=(-((-(y(1)*(1-params(1))))/(y(4)*y(4))));
  g1(4,7)=1;
  g1(5,1)=1;
  g1(5,3)=(-(T48*y(5)*getPowerDeriv(y(3),params(1),1)));
  g1(5,4)=(-(T47*getPowerDeriv(y(4),1-params(1),1)));
  g1(5,5)=(-(T46*T48));
  g1(6,3)=1-(1-params(3));
  g1(6,8)=(-1);
  g1(7,1)=1;
  g1(7,2)=(-1);
  g1(7,8)=(-1);
  g1(8,5)=1/y(5)-params(6)*1/y(5);
  g1(9,9)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,81);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,729);
end
end
end
end
