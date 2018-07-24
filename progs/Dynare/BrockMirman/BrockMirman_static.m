function [residual, g1, g2, g3] = BrockMirman_static(y, x, params)
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

residual = zeros( 3, 1);

%
% Model equations
%

T9 = y(1)^(-1);
T18 = y(2)^(params(1)-1);
T21 = y(2)^params(1);
lhs =T9;
rhs =T9*params(1)*params(2)*y(3)*T18;
residual(1)= lhs-rhs;
lhs =y(2);
rhs =y(3)*T21-y(1);
residual(2)= lhs-rhs;
lhs =log(y(3));
rhs =log(y(3))*params(3)+params(4)*x(1);
residual(3)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(3, 3);

  %
  % Jacobian matrix
  %

T33 = getPowerDeriv(y(1),(-1),1);
  g1(1,1)=T33-T18*y(3)*params(1)*params(2)*T33;
  g1(1,2)=(-(T9*params(1)*params(2)*y(3)*getPowerDeriv(y(2),params(1)-1,1)));
  g1(1,3)=(-(T9*params(1)*params(2)*T18));
  g1(2,1)=1;
  g1(2,2)=1-y(3)*getPowerDeriv(y(2),params(1),1);
  g1(2,3)=(-T21);
  g1(3,3)=1/y(3)-params(3)*1/y(3);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],3,9);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],3,27);
end
end
end
end
