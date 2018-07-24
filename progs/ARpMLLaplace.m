function ML = ARpMLLaplace(y,p,const,alph)
% =======================================================================
% Maximum Likelihood Estimation of Laplace AR(p) model:
% y_t = c + d*t + theta_1*y_{t-1} + ... + theta_p*y_{t-p} + u_t
% with u_t ~ Laplace distributed with E(u_t)=0, Var(u_t)=2
% =======================================================================
% ML = ARpMLLaplace(y,p,const,alph)
% -----------------------------------------------------------------------
% INPUTS
%   - y     : dependent variable vector. [Tx1]
%   - p     : number of lages. [scalar]
%   - const : 0 no constant; 1 constant; 2 constant and linear trend. [scalar]
%   - alph  : significance level for t statistic. [scalar]
% -----------------------------------------------------------------------
% OUTPUT
%	- ML: structure including estimation results
%     - T_eff          : effective sample size used in estimation. [scalar]
%     - thetatilde     : estimate of coefficients. [(const+p)x1]
%     - sig_thetatilde : estimate of standard error of coefficients. [(const+p)x1]
%     - tstat          : t statistics given alph as significance level. [(const+p)x1]
%     - pvalues        : p values of H_0: thetahat = 0. [(const+p)x1]
%     - sig_utilde     : estimate of standard deviation of error term u. [scalar]
%     - theta_ci       : (1-alph)% confidence intervall for theta given
%                        significance level alph. [(const+p)x2]
%     - logl           : value of maximized log likelihood. [scalar]
% -----------------------------------------------------------------------
% CALLS
%    - LogLikeARpLaPlace : Computes log likelihood function of Laplace AR(p)
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================

T = size(y,1);                      % sample size
x0 = randn(p+const,1);              % randomize start values, note that sig_u is known

% Optimization with fminunc which finds the minimum of negative log-likelihood
f = @(x)-1*LogLikeARpLaplace(x,y,p,const); % use function handle to hand over 
                                           % additional parameters to negative
                                           % of LogLikeARpLaplace
[x,fval,exitflag,output,grad,hessian] = fminunc(f,x0);

thetatilde = x;                     % estimated coefficient values
V = inv(hessian);                   % estimated covariance matrix of coefficients
sig_thetatilde = sqrt(diag(V));     % estimated standard error of coefficients
T_eff = T-p;                        % effective sample size used in estimation
logl  = -fval;                      % value of maximized log likelihood
tstat = thetatilde./sig_thetatilde; % t-statistic
tcrit = -tinv(alph/2,T_eff-p);      % critical value from t-distribution
pvalues = tpdf(tstat,T_eff-p);      % p-values from t-distribution

% confidence interval given signficance level alph
theta_ci=[thetatilde-tcrit.*sig_thetatilde, thetatilde+tcrit.*sig_thetatilde];

% Store into output structure
ML.T_eff          = T_eff;
ML.logl           = logl;
ML.thetatilde     = thetatilde;
ML.sig_thetatilde = sig_thetatilde;
ML.tstat          = tstat;
ML.pvalues        = pvalues;
ML.theta_ci       = theta_ci;

end %function end