%   VAR Estimation by OLS Method with Lagged Exogenous Variables
%   FUNCTION : 
%   [Psi, Alpha, PhiL, R, Res, StErr, Sigma] = VARs(y, x, ylag, xlag, deter)
% 
%   INPUT :  i) DY    :   Endogenous Variables      (T, m) 
%           ii) ECM   :   Eqquilibrium Relationship (T, h) > Y(T-1) * Beta
%          iii) ylag  :   Number of lags for DY
%           iv) deter :   Deterministic Variables   (T, ndet)
%
%   OUTPUT :  i) Psi ; Alpha  : Constant and Exogenous Coefs
%            ii) PhiL         : Endog. Lagged Coefficients (neqs, neqs*ylag )
%           iii) Res          : Residuals
%             v) SEBeta       : Standard Errors of Regressors
%            vi) Sigma        : Residual Covariance Matrix (Unbiased): SSR/(AdjNObs - # Parameters)
%
% Agostino © Oct 2005

function [Psi, Alpha, PhiL, Res, SEBeta, Sigma] = VECM(dy, ecm, ylag, deter)
   
    % >> LAGGING the Vector y and x
        [T, m ]   = size(dy);
        [Ylag, Y] = LagByLags(dy, ylag);
        X         = ecm( ylag + 1 : end, : );
        nexo      = size(X, 2);

	% >> DETERMINISTIC COMPONENT 
        if     and( eq(length(deter), 1), eq(deter, 1) )
            deter = ones( T - ylag , 1);
        elseif and( eq(length(deter), 1), eq(deter, 0) )
            deter = [];
        else
            deter = deter(ylag + 1 : end, :);
        end
            ndet  = size(deter, 2);
    
    % >> ORDERING : DETER | ECM | Endogenous
        R.ORDER  = ' DETER - X - Ylag';
        WeakExog = [ deter, X, Ylag ];
        epar     = size(WeakExog, 2);

    % >> ESTIMATION & INFERENCE    
        Coefs = WeakExog\Y;
        Res   = Y - WeakExog * Coefs;
        Sigma = ( Res'*Res ) / ( T - ylag - epar );    % DoF Adjusted * !!! This formula is correct !!! *
        XXI   = (WeakExog'*WeakExog)\eye(epar);
        VBeta = mat( diag( kron( Sigma , XXI ) ), m );
        SEBeta= sqrt(VBeta);
        TStat = Coefs./SEBeta; 
        
    % >> PRINTING MATRICES
        Psi   = Coefs(1:ndet, :)';
        Alpha = Coefs(ndet + 1 : ndet + nexo, :)';
        PhiL  = Coefs(ndet + nexo + 1 : end , :)';
        if not( isempty(PhiL))
            Phi3D = Pages(PhiL);
        else
            PhiL  = [];
            Phi3D = [];
        end
        