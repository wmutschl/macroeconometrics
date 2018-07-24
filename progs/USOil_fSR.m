function f=USOil_fSR(B0inv,SIGMAUHAT,quickanddirty)
% =======================================================================
% Evaluates the system of nonlinear equations vech(SIGMAUHAT) = vech(B0inv*B0inv')
% subject to the short-run restrictions.
% =======================================================================
% f=USOil_fSR(B0inv,SIGMAUHAT,quickanddirty)
% -----------------------------------------------------------------------
% INPUTS
%   - B0inv     : candidate for short-run impact matrix. [nvars x nvars]
%   - SIGMAUHAT : covariance matrix of reduced-form residuals. [nvars x nvars]
%   - quickanddirty: optional flag to define restrictions in a more graphical way, see below
% -----------------------------------------------------------------------
% OUTPUTS
%   - f : function value, see below
% -----------------------------------------------------------------------
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu
if nargin < 3
    quickanddirty = 0;
end

if quickanddirty == 1
    f=[vech(B0inv*B0inv'-SIGMAUHAT);
        B0inv(1,2) - 0;
        B0inv(1,3) - 0;
        B0inv(2,3) - 0;
        ];
else
    % nan means unconstrained, 0 (or any other number) restricts to this number
    Rshort = [nan   0   0;
              nan nan   0;
              nan nan nan;
              ];
    selSR  = find(isnan(Rshort)==0); % Index for short-run restrictions on impact
    f=[vech(B0inv*B0inv'-SIGMAUHAT);
       B0inv(selSR) - Rshort(selSR);
       ];
end