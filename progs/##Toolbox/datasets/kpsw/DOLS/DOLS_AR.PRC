
proc(1) = dols_ar(y,x,fperg,lperg,nlag,nar);

/*  dols_ar.prc, gmg 11/06/90, JHS 11/7/90
    This procedure calculates cointegrating vectors by dynamic gls.
    Can be used on any system of equations or a single equation.

 Input:

y = y data (a matrix of dep. variables)
x = x data (a matrix of regressors, in levels)
fperg = first obs of y and x to use in regression
        (must have available at least nlag+nar+1 observations before fperg)
lperg = last obs of y and x to use in regression
        (must have available at least nlag observations after lperg)
nlag = # of lags
nar = # of AR lags
*/

local *;
k = cols(y);         @ Number of equations in VAR @
nv = cols(x);        @ Number of variables in the right hand side (levels) @
nrhsv = 1+2*nv+2*nlag*nv; @ Number of right hand side variables (inc
                            constant, leads and lags)@

xa=x[fperg:lperg,.];
ya=y[fperg:lperg,.];
" ****** RESULTS of DOLS system estimation ******";" ";

/* Add constant and lagged differences to X matrix */

i=nlag; do while i>= -nlag;
xa=xa~(x[fperg+i:lperg+i,.]-x[fperg+i-1:lperg+i-1,.]);
i=i-1;endo;
xa=ones(rows(xa),1)~xa; @ filtering, also, the constant term @

@ DOLS regression:  also serves as initial parameter vector for DGLS  @

xxi=invpd(xa'xa);
xy=xa'ya;
beta=xxi*xy;

@ DOLS regression: sample from fperg to lperg @

@xa2=xa[nar+1:rows(xa),.];
ya2=ya[nar+1:rows(ya),.];
xxi2=invpd(xa2'xa2);
xy=xa2'ya2;
beta2=xxi2*xy;
@
vecb=vec(beta[2:(nv+1),.]);

@ calculate phihat @
kr=ya-xa*beta;
xr=kr[nar:rows(kr)-1,.];
i=2; do while i <= nar;
  xr=xr~kr[nar+1-i:rows(kr)-i,.];
i=i+1;endo;
rho=(xr'kr[nar+1:rows(kr),.])/(xr'xr);
phihat=rho';
res=kr[(nar+1):rows(kr),.]-xr*rho;
sigma=(res'res)/(rows(res)-nrhsv);

@ form phi(1) matrix @
phione=eye(k);
i=1; j1=1; j2=k; do while i <= nar;
  phione=phione-phihat[.,j1:j2]; j1=j1+k; j2=j2+k;
i=i+1; endo;

@ Calculate variance of bet @
spec=inv(phione)*sigma*(inv(phione))';
varb=spec.*.xxi;

@ form relevant varb and se's @
indx=reshape(seqa(1,1,k*nrhsv),k,nrhsv);indx=indx';
indx=vec(indx[2:(nv+1),.]);
relvb=varb[indx,indx];
relse=sqrt(diag(relvb));

@ ------ output results ---------- @
format /RDN 10,0;
" "; " --------- Dynamic OLS estimates of cointegrating vector ---------";
"    first, last period of regression:" fperg lperg;
"        no. of leads (= no. of lags):" nlag;
" order of Cochrane-Orcutt correction:" nar;
format /RDN 11,5;" ";
 "   Matrix of point estimates:"; reshape(vecb,cols(y),cols(x));" ";
 "   Matrix of standard errors:"; reshape(relse,cols(y),cols(x));" ";
 " PHIHAT: VAR lag polynomial matrix for errors in DOLS equations,";
 "         ordered equation by equation:";
 phihat;

retp(vecb~relvb);
endp;
