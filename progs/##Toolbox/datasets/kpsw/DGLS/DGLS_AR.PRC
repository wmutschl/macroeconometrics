
proc(1) = dgls_ar(y,x,fperg,lperg,nlag,nar,maxit);

/*  dgls_ar.prc, gmg 10/31/90
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
maxit = # of iterations
*/

local *;
k = cols(y);         @ Number of equations in VAR @
nv = cols(x);        @ Number of variables in the right hand side (levels) @
nrhsv = 1+2*nv+2*nlag*nv; @ Number of right hand side variables (inc
                            constant, leads and lags)@

xa=x[fperg-nar:lperg,.];
ya=y[fperg-nar:lperg,.];

/* Add constant and lagged differences to X matrix */

i=nlag; do while i>= -nlag;
xa=xa~(x[fperg-nar+i:lperg+i,.]-x[fperg-nar+i-1:lperg+i-1,.]);
i=i-1;endo;
xa=ones(rows(xa),1)~xa; @ filtering, also, the constant term @

@ OLS  @

xxi=invpd(xa'xa);
xy=xa'ya;
beta=xxi*xy;
" ****** RESULTS of DGLS system estimation ******";" ";
format /RDN 10,6;
" DOLS estimate of cointegrating vector: ";
beta[2:(nv+1),.]';

@ corc iterations @

it=1; do while it <= maxit;
kr=ya-xa*beta;
xr=kr[nar:rows(kr)-1,.];
i=2; do while i <= nar;
xr=xr~kr[nar+1-i:rows(kr)-i,.];
i=i+1;endo;
rho=(xr'kr[nar+1:rows(kr),.])/(xr'xr);
phihat=rho';

@ filter all variables @

yar=ya[nar+1:rows(ya),.];
xar=eye(k).*.xa[nar+1:rows(xa),.];
i=1; j1=1; j2=k;
do while i <= nar;
 yar=yar-ya[nar+1-i:rows(ya)-i,.]*rho[j1:j2,.];
 xar=xar-phihat[.,j1:j2].*.xa[nar+1-i:rows(xa)-i,.];
 j1=j1+k;
 j2=j2+k;
 i=i+1;
endo;

@ form ytil @
ytil=yar[.,1];
i=2; do while i <= k;
ytil=ytil|yar[.,i];
i=i+1; endo;

@ regress ytil on xar @
xxi=invpd(xar'xar);
xy=xar'ytil;
bet=xxi*xy;
@beta=reshape(bet,nrhsv,k);beta=beta';@
beta=bet[1:nrhsv,.];
j=nrhsv; i=2;
do while i <= k;
beta=beta~bet[j+1:j+nrhsv,.];
j=j+nrhsv;i=i+1;
endo;
it=it+1; endo;

s2=(ytil'ytil-bet'(xar'ytil))/(rows(ytil)-k*nrhsv);
varb=s2*xxi;

@ form vec(relevant beta's) and relevant varb and se's @
/*
vecb=1;indx=1;
j2=nv+1; j=2; i=1;
do while i <= k;
  do while j <= j2;
     indx=indx~j;
     vecb=vecb|bet[j,.];
     j=j+1;
   endo;
j=2+i*nrhsv; j2=nv+1+i*nrhsv; i=i+1;
endo;
vecb=vecb[2:rows(vecb),.];
indx=indx[.,2:cols(indx)]; */

vecb=vec(beta[2:nv+1,.]);
indx=reshape(seqa(1,1,k*nrhsv),k,nrhsv);indx=indx';
indx=vec(indx[2:(nv+1),.]);
relvb=varb[indx,indx];
relse=sqrt(diag(relvb));
" ";"s2= ";; s2; " ";

@ ------ output results ---------- @
format /RDN 10,6;
" PHIHAT: VAR lag polynomial matrix for errors in DGLS equations,";
"         ordered equation by equation -- final estimates";
"         in Cochrane-Orcutt correction:";
phihat;

" "; " --------- Dynamic GLS estimates of cointegrating vector ---------";
 "     vecb'"; vecb';
 "   std errors"; relse';" ";

retp(vecb~relvb);
endp;
