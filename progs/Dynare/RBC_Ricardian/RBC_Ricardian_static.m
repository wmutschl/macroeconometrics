function [residual, g1, g2, g3] = RBC_Ricardian_static(y, x, params)
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

residual = zeros( 14, 1);

%
% Model equations
%

T58 = (y(7)/y(9))^params(1);
T63 = (y(9)/y(7))^(1-params(1));
T67 = y(7)^params(1);
T68 = y(14)*T67;
T69 = y(9)^(1-params(1));
lhs =1;
rhs =params(2)*(1-params(3)+y(13));
residual(1)= lhs-rhs;
lhs =y(4);
rhs =y(12)*y(11);
residual(2)= lhs-rhs;
lhs =y(3);
rhs =y(12)*params(4)/(1-params(4))*(1-y(10));
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(12)*params(4)/(1-params(4))*(1-y(11));
residual(4)= lhs-rhs;
lhs =y(9);
rhs =y(10)*params(5)+y(11)*(1-params(5));
residual(5)= lhs-rhs;
lhs =y(2);
rhs =y(3)*params(5)+y(4)*(1-params(5));
residual(6)= lhs-rhs;
lhs =y(7);
rhs =params(5)*y(8);
residual(7)= lhs-rhs;
lhs =y(5);
rhs =params(5)*y(6);
residual(8)= lhs-rhs;
lhs =y(12);
rhs =(1-params(1))*y(14)*T58;
residual(9)= lhs-rhs;
lhs =y(13);
rhs =params(1)*y(14)*T63;
residual(10)= lhs-rhs;
lhs =y(1);
rhs =T68*T69;
residual(11)= lhs-rhs;
lhs =y(8);
rhs =y(6)+(1-params(3))*y(8);
residual(12)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(5);
residual(13)= lhs-rhs;
lhs =log(y(14));
rhs =log(y(14))*params(6)+x(1);
residual(14)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(14, 14);

  %
  % Jacobian matrix
  %

T86 = getPowerDeriv(y(7)/y(9),params(1),1);
T93 = getPowerDeriv(y(9)/y(7),1-params(1),1);
  g1(1,13)=(-params(2));
  g1(2,4)=1;
  g1(2,11)=(-y(12));
  g1(2,12)=(-y(11));
  g1(3,3)=1;
  g1(3,10)=(-(y(12)*(-(params(4)/(1-params(4))))));
  g1(3,12)=(-(params(4)/(1-params(4))*(1-y(10))));
  g1(4,4)=1;
  g1(4,11)=(-(y(12)*(-(params(4)/(1-params(4))))));
  g1(4,12)=(-(params(4)/(1-params(4))*(1-y(11))));
  g1(5,9)=1;
  g1(5,10)=(-params(5));
  g1(5,11)=(-(1-params(5)));
  g1(6,2)=1;
  g1(6,3)=(-params(5));
  g1(6,4)=(-(1-params(5)));
  g1(7,7)=1;
  g1(7,8)=(-params(5));
  g1(8,5)=1;
  g1(8,6)=(-params(5));
  g1(9,7)=(-((1-params(1))*y(14)*1/y(9)*T86));
  g1(9,9)=(-((1-params(1))*y(14)*T86*(-y(7))/(y(9)*y(9))));
  g1(9,12)=1;
  g1(9,14)=(-((1-params(1))*T58));
  g1(10,7)=(-(params(1)*y(14)*(-y(9))/(y(7)*y(7))*T93));
  g1(10,9)=(-(params(1)*y(14)*T93*1/y(7)));
  g1(10,13)=1;
  g1(10,14)=(-(params(1)*T63));
  g1(11,1)=1;
  g1(11,7)=(-(T69*y(14)*getPowerDeriv(y(7),params(1),1)));
  g1(11,9)=(-(T68*getPowerDeriv(y(9),1-params(1),1)));
  g1(11,14)=(-(T67*T69));
  g1(12,6)=(-1);
  g1(12,8)=1-(1-params(3));
  g1(13,1)=1;
  g1(13,2)=(-1);
  g1(13,5)=(-1);
  g1(14,14)=1/y(14)-params(6)*1/y(14);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,2744);
end
end
end
end
