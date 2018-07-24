function f=RWZSRLR_f(B0inv,SIGMAUHAT,LRMat)
% =======================================================================
% Evaluates the system of nonlinear equations vech(SIGMAUHAT) = vech(B0inv*B0inv')
% subject to the short-run and long-run restrictions.
% =======================================================================
% f=RWZSRLR_f(B0inv,SIGMAUHAT,LRMat)
% -----------------------------------------------------------------------
% INPUTS
%   - B0inv     : candidate for short-run impact matrix. [nvars x nvars]
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
%	- LRMat     : total long-run impact matrix. [nvars x nvars] 
%                 - for VAR model A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value, see below
% -----------------------------------------------------------------------
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu

LRMAT = LRMat*B0inv;

f=[vech(SIGMAUHAT-B0inv*B0inv');
        B0inv(1,1) - 0;
        LRMAT(1,1) - 0;
        LRMAT(1,2) - 0;
        ];

end