function loglik=LogLikeARpNorm(theta_0,y,p,const)
%likelihood function to estimate AR(p) model with LaPlace distributed errors
theta0 = theta_0(1:(const+p));
[T, K] = size(y);
T_eff = T-p;
Y = lagmatrix(y,1:p);
if const==1 %constant
    Y = [ones(T_eff,1) Y((p+1):end,:)];
elseif const==2 % time trend and constant
    Y = [ones(T_eff,1) 1:T_eff' Y((p+1):end,:)];
else
    Y = Y((p+1):end,:);
end
y = y(p+1:end);
uhat = y - Y*theta0;

loglik = -T_eff*log(2)-sum(abs(uhat));
 
if isnan(loglik) || isinf(loglik) || ~isreal(loglik)
    loglik = -1e10;
end

end % function end