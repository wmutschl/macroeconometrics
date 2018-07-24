function f=gali1992_f(B0inv,SIGMAUHAT,LRMat)
% =======================================================================
% Evaluates the system of nonlinear equations vech(SIGMAUHAT) = vech(B0inv*B0inv')
% subject to the short-run and long-run restrictions.
% =======================================================================
% f=RWZSRLR_f(B0inv,SIGMAUHAT,LRMat)
% -----------------------------------------------------------------------
% INPUTS
%   - B0_diageps : [nvars x (nvars+1)] candidate matrix for short-run impact matrix [nvars x nvars] and diagonal elements of SIGeps (nvar x 1) 
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
%	- LRMat     : total long-run impact matrix. [nvars x nvars] 
%                 - for VAR model A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value, see below
% -----------------------------------------------------------------------
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
%B0inv = inv(B0); % Get short-run impact matrix
LRMAT = LRMat*B0inv;

f=[vech(SIGMAUHAT-B0inv*B0inv');
   B0inv(1,2) - 0;
   B0inv(1,3) - 0;
   B0inv(2,3) + B0inv(2,4) - 0;
   LRMAT(1,2) - 0;
   LRMAT(1,3) - 0;
   LRMAT(1,4) - 0;
   ];
% f=[vech(B0*SIGMAUHAT*B0'-eye(size(B0)));
%    B0(1,2) - 0;
%    B0(1,3) - 0;
%    B0(2,3) + B0(2,4) - 0;
%    %B(2,1) - 0;
%    %B0(3,3) - 0;
%    LRMAT(1,2) - 0;
%    LRMAT(1,3) - 0;
%    LRMAT(1,4) - 0;
%    ];
end