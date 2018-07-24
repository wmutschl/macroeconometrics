function [residual, g1, g2, g3] = RBC_habit_static(y, x, params)
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

cback__ = 1/(y(2)-y(2)*params(6));
c__ = 1/(y(2)-y(2)*params(6));
cp__ = y(9);
T49 = (y(4)/y(5))^params(1);
T54 = (y(5)/y(4))^(1-params(1));
T58 = y(4)^params(1);
T59 = y(8)*T58;
T60 = y(5)^(1-params(1));
lhs =(cback__-params(6)*params(2)*c__)/(c__-params(6)*params(2)*cp__);
rhs =params(2)*(1+y(7)-params(3));
residual(1)= lhs-rhs;
lhs =(cback__*params(4)-c__*params(6)*params(2)*params(4))*y(6);
rhs =(1-params(4))/(1-y(5));
residual(2)= lhs-rhs;
lhs =y(6);
rhs =(1-params(1))*y(8)*T49;
residual(3)= lhs-rhs;
lhs =y(7);
rhs =params(1)*y(8)*T54;
residual(4)= lhs-rhs;
lhs =y(1);
rhs =T59*T60;
residual(5)= lhs-rhs;
lhs =y(4);
rhs =y(4)*(1-params(3))+y(3);
residual(6)= lhs-rhs;
lhs =y(1);
rhs =y(2)+y(3);
residual(7)= lhs-rhs;
lhs =log(y(8));
rhs =log(y(8))*params(5)+x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =1/(y(2)-y(2)*params(6));
residual(9)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

  %
  % Jacobian matrix
  %

T80 = (-(1-params(6)))/((y(2)-y(2)*params(6))*(y(2)-y(2)*params(6)));
T94 = getPowerDeriv(y(4)/y(5),params(1),1);
T101 = getPowerDeriv(y(5)/y(4),1-params(1),1);
  g1(1,2)=((c__-params(6)*params(2)*cp__)*(T80-params(6)*params(2)*T80)-(cback__-params(6)*params(2)*c__)*T80)/((c__-params(6)*params(2)*cp__)*(c__-params(6)*params(2)*cp__));
  g1(1,7)=(-params(2));
  g1(1,9)=(-((cback__-params(6)*params(2)*c__)*(-(params(6)*params(2)))))/((c__-params(6)*params(2)*cp__)*(c__-params(6)*params(2)*cp__));
  g1(2,2)=y(6)*(params(4)*T80-params(6)*params(2)*params(4)*T80);
  g1(2,5)=(-((1-params(4))/((1-y(5))*(1-y(5)))));
  g1(2,6)=cback__*params(4)-c__*params(6)*params(2)*params(4);
  g1(3,4)=(-((1-params(1))*y(8)*1/y(5)*T94));
  g1(3,5)=(-((1-params(1))*y(8)*T94*(-y(4))/(y(5)*y(5))));
  g1(3,6)=1;
  g1(3,8)=(-((1-params(1))*T49));
  g1(4,4)=(-(params(1)*y(8)*(-y(5))/(y(4)*y(4))*T101));
  g1(4,5)=(-(params(1)*y(8)*T101*1/y(4)));
  g1(4,7)=1;
  g1(4,8)=(-(params(1)*T54));
  g1(5,1)=1;
  g1(5,4)=(-(T60*y(8)*getPowerDeriv(y(4),params(1),1)));
  g1(5,5)=(-(T59*getPowerDeriv(y(5),1-params(1),1)));
  g1(5,8)=(-(T58*T60));
  g1(6,3)=(-1);
  g1(6,4)=1-(1-params(3));
  g1(7,1)=1;
  g1(7,2)=(-1);
  g1(7,3)=(-1);
  g1(8,8)=1/y(8)-params(5)*1/y(8);
  g1(9,2)=(-T80);
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
