function [vecy]=vec(y)

% Lutz Kilian
% University of Pennsylvania
% December 1993

% This function vectorizes an (a x b) matrix y.  The resulting vector vecy
% has dimension (a*b x 1).

[row,column]=size(y);
vecy=reshape(y,row*column,1);
