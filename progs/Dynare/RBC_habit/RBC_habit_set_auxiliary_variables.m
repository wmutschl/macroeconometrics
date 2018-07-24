function y = RBC_habit_set_auxiliary_variables(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

y(9)=1/(y(2)-y(2)*params(6));
