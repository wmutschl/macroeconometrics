function ML = ARpML(y,p,const,distrib)
alph = 0.5;
[T, K] = size(y);
if distrib == 0
    x0 = randn(p+const+1,1); % randomize start values
    f = @(x)-1*LogLikeARpNorm(x,y,p,const); % use function handle to handle additional parameters to LogLikeARpNorm
elseif distrib == 1
    x0 = randn(p+const,1); % randomize start values
    f = @(x)-1*LogLikeARpLaPlace(x,y,p,const); % use function handle to handle additional parameters to LogLikeARpLaPlace
end
[x,fval,exitflag,output,grad,hessian] = fminunc(f,x0);
thetatilde = x(1:p+const);
V = inv(hessian);
sig = sqrt(diag(V));
sig_thetatilde = sig(1:p+const);
T_eff = T-p;
logl = -fval;

tstat = thetatilde./sig_thetatilde;
tcrit=-tinv(alph/2,T_eff-p*K);
pvalues = tpdf(tstat,T_eff-p*K);
theta_ci=[thetatilde-tcrit.*sig_thetatilde, thetatilde+tcrit.*sig_thetatilde];

ML.T_eff = T_eff;
ML.logl = logl;
ML.thetatilde = thetatilde;
ML.sig_thetatilde = sig_thetatilde;
ML.tstat = tstat;
ML.pvalues = pvalues;
ML.theta_ci = theta_ci;
end %function end