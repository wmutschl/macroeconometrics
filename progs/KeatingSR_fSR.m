function f=KeatingSR_fSR(B0_diageps,SIGMAUHAT)
% =======================================================================
% Evaluates the system of nonlinear equations vech(SIGMAUHAT) =
% vech(B0inv*SIGeps*B0inv') where SIGeps has only values on the diagonal
% given a candidate for B0 and the diagonal elements of SIGeps
% subject to the short-run restrictions in the Keating (1992) model.
% =======================================================================
% f=keatingSR_f_SR(B0_diageps,SIGMAUHAT)
% -----------------------------------------------------------------------
% INPUTS
%   - B0_diageps : [nvars x (nvars+1)] candidate matrix for short-run impact matrix [nvars x nvars] and diagonal elements of SIGeps (nvar x 1) 
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value, see below
% -----------------------------------------------------------------------
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu

nvars = size(SIGMAUHAT,1);  % Number of variables
B0 = B0_diageps(:,1:nvars); % Get B0 matrix from candidate
diageps = diag(B0_diageps(1:nvars,nvars+1)); % Get diagonal elements from candidate, make diagonal matrix
B0inv = inv(B0); % Get short-run impact matrix

f=[vech(B0inv*diageps*B0inv'-SIGMAUHAT);
   B0(1,1) - 1;
   B0(3,1) - 0;
   B0(4,1) - B0(4,2);
   B0(1,2) - 0;
   B0(2,2) - 1;
   B0(3,2) - 0;
   B0(1,3) - 0;
   B0(3,3) - 1;
   B0(1,4) - 0;
   B0(4,4) - 1;
   ];