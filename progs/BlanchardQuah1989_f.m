function f=BlanchardQuah1989_f(B0inv,SIGMAUHAT,LRMat)
% =======================================================================
% Evaluates the system of nonlinear equations 
% vech(SIGMAUHAT) = vech(B0inv*B0inv')
% subject to the long-run restrictions
% [* 0;
%  * *];
% =======================================================================
% f=BlanchardQuah1989_f(B0inv,SIGMAUHAT,LRMat)
% -----------------------------------------------------------------------
% INPUTS
%   - B0inv     : candidate for short-run impact matrix. [nvars x nvars]
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
%	- LRMat     : total long-run impact matrix. [nvars x nvars] 
%                 - for VAR model A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value of system of nonlinear equations
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu
% =======================================================================

THETA = LRMat*B0inv; % Cumulated (long-run) impulse response function
f=[vech(B0inv*B0inv'-SIGMAUHAT);
    THETA(1,2) - 0;
    ];
end