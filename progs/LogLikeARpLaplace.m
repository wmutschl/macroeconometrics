function loglik=LogLikeARpLaplace(x,y,p,const)
% =======================================================================
% Computes the conditional log likelihood function of Laplace AR(p) model:
% y_t = c + d*t + theta_1*y_{t-1} + ... + theta_p*y_{t-p} + u_t
% with u_t ~ Laplace distributed with E(u_t)=0, Var(u_t)=2
% =======================================================================
% loglik=LogLikeARpLaplace(x,y,p,const)
% -----------------------------------------------------------------------
% INPUT
%   - x     : [c,d,theta_1,...,theta_p]'. [(const+p)x1]
%   - y     : data vector of dimension T. [Tx1]
%   - p     : number of lags. [scalar]
%   - const : 1 constant; 2 constant and linear trend in model. [scalar]
% -----------------------------------------------------------------------
% OUTPUT
%	- loglik : value of Laplace log-likelihood
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================

theta = x; % AR coefficients
T = size(y,1);          % sample size

Y = lagmatrix(y,1:p);   % create matrix with lagged variables
if const == 1           % add constant
    Y = [ones(T,1) Y];
elseif const == 2       % add constant and time trend
    Y = [ones(T,1) transpose(1:T) Y];
end
Y = Y((p+1):end,:);     % get rid of initial observations
y = y(p+1:end);         % get rid of initial observations

uhat = y - Y*theta;     % ML residuals

% compute the conditional log likelihood
loglik = -log(2)*(T-p) -sum(abs(uhat));

if isnan(loglik) || isinf(loglik) || ~isreal(loglik)
    loglik = -1e10;     % if anything goes wrong set value very small
end

end % function end