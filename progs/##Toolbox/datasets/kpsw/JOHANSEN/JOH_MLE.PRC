@========================== JOH_MLE.PRC ==============================
   modification of KPSW1.prg to include constant, time trend correctly
   This file includes main proc, JOH_MLE, and canonical correlation
   proc, CANON @


@ Load test statistics table. @
@chkf= files(sfname,0); if chkf[1,1] $== 0; goto stop2; endif;
stop2: fname= sfname; goto emsg;
emsg: "ERROR: THE FILE " upper(fname) " CANNOT BE FOUND.";@
load joh_n[399,7]= "\\kpsw\\johansen\\joh.tab";
load joh_c[399,7]= "\\kpsw\\johansen\\joh_c.tab";
load joh_t[399,7]= "\\kpsw\\johansen\\joh_t.tab";

@======================================================================
Procedure to calculate canonical vectors and correlations,
based on the algorithm in Kennedy & Gentle(1980), Ch.12.2, p.565
Also see Anderson(1984), Ch.12.3.2, p.495-7.

Canonical variates, or cointegrating vector, "beta",
are return variables by this procedure.
squared canonical correlations, "lmbd", will be returned
as a global variable.
======================================================================@
proc(2)=canon(rk,r0,r,p);
@local zk, z0, zkk, z00, omga, beta, mb;
clearg lmbd;@
local *;

@ form the adjusted data matrices zk & z0 @
@" **** means of residuals:"; meanc(rk)';meanc(r0)';@

zk= rk - meanc(rk)';  zkk= zk'zk./rows(zk);
z0= r0 - meanc(r0)';  z00= z0'z0./rows(z0);

@ factor zk & z0 by QR decomposition @
clearg zkq, zkr, z0q, z0r;

{zkq,zkr}=qr(zk);
{z0q,z0r}=qr(z0);

@ singular-value decompositon of zkq'z0q @
clearg svdu, svdv;

omga= zkq'z0q;
{svdu,lmbd,svdv}=svd1(omga);     @ lmbd is a canon. corr. @
lmbd= diag(lmbd).^2;             @ eigen values are squared canon. corr. @

beta= inv(zkr)*svdu;

@ normalize beta such that beta'*zkk*beta=1 @
mb= sqrt( 1./diag(beta'zkk*beta) );
beta= beta.*mb';

@-----------------------------
 "... optimal beta can be found from beta(hat)
   by beta=beta(hat)*rho for rho an r*r matrix of full rank"
 Johansen(1987), p.8
------------------------------@
beta= beta*diagrv( zeros(cols(beta),cols(beta)),1/beta[1,.] );

retp(beta,lmbd);
endp;
@::::::::::::::::::::::: End of procedure ::::::::::::::::::::::@


@========================== JOH_MLE.PRC ==============================@
proc(1) = joh_mle(y,x,fperg,lperg,nlag,optn,ah,ep,ictst);

/*
 Input:

y = y data (a matrix of dep. variables)
x = x data (a matrix of regressors, in levels)
     Note: cols(y) = number of cointegrating vectors
fperg = first obs of y and x to use in regression
        (must have available at least nlag+nar+1 observations before fperg)
lperg = last obs of y and x to use in regression
        (must have available at least nlag observations after lperg)
nlag = # of lags
nar = # of AR lags
optn = 1: no deterministic terms in VECM regressions
       2: constant in VECM regressions
       3: constant, time in VECM regressions
*/
 @ loadp qr, svd, timestr, datestr;@
 local *;
 " ";
"=====================================================================";
"---- Johansen estimator of cointegrating vectors, test statistcs ----";
"=====================================================================";
format /RDN 5,0;"first, last period of regression:" fperg lperg;
format /RDN 11,5;

xx=y~x;xx=xx[(fperg-nlag-1):lperg,.];
r=cols(y); p=cols(xx);

k=nlag+1;
@::::: construct VAR system of equation (2.3) in p.6 :::::@
 xlg= xx[2:rows(xx),.]-xx[1:(rows(xx)-1),.];
 lfx= xlg[k:(rows(xx)-1),.];                  @ dX(t) @
 xk= xx[1:(rows(xx)-k),.];                    @ lagged level: X(t-k) @
 rtx= zeros(rows(xk),1 );
 nperg=rows(xk);

 t=1; do while t<=(k-1);
  rtx= rtx~xlg[t:(rows(xx)-k+t-1),.];
 t=t+1; endo;


 rtx= rtx[.,2:cols(rtx)];               @ lagged differences @
 "Treatment of deterministic components = ";;
 if optn==1;"None";tab=joh_n;
 elseif optn==2;"Constant in Regression";tab=joh_c;
  rtx=ones(rows(rtx),1)~rtx;
 elseif optn==3;"Constant, time in Regression";tab=joh_t;
  rtx=ones(rows(rtx),1)~seqa(1,1,rows(rtx))~rtx;
 endif;

@----------------------------------------
regress X(t-k) & dX(t) on the lagged differences
  and get residuals R(kt) & R(0t)

S(i,j) in Johansens's will be ( R(i)'R(j) )./( rows(R) )
----------------------------------------@

xxm= moment( rtx, 0 );
rk= xk - rtx*solpd( rtx'xk, xxm );
r0= lfx - rtx*solpd( rtx'lfx, xxm );

@----------------------------------------
Johansen's Theorem I:
 Test statistics that there are at most r cointegrating vectors
----------------------------------------@
ir=1;
do until ir>2;
r=(ir==2)*cols(y);
{beta,lmbd1} = canon(rk,r0,r,p);

if r .>= 1;
 @ choose first r canonical vectors @
 beta= beta[.,1:r];
 " ";"  Cointegrating vector, Johansen normalization:"; beta;

   @ Bewley Form @

 inm=rows(beta)-cols(beta);
@ wat=beta[inm+1:rows(beta),.];@
 wat=beta[1:cols(beta),.];
 beta=beta*inv(wat);
format /RDN 11,5;
 " ";"  Cointegrating vector -- Bewley Form:"; beta;
endif;

@ test statistics @
" ";"canonical correlations";lmbd1;
tst1 = -nperg*sumc( ln( 1-lmbd1[r+1:p,1] ) );
df = r; pv= 100-( counts(tab[.,p-r+1],tst1) )./4;

" ";"Test statistic that there are at most ";; format /RDN 2,0; r;;
    " cointegrating vectors: "
" ";"   Test statistic: ";;format /RDN 7,3; tst1;;
"       D.F.= ";;format /RDN 3,0;df;;
"       P-value(%)= ";; format /RDN 4,2; pv;format /RDN 11,5;
ir=ir+1; endo;

/* carry out test of cointegrating vectors if ictst = 1 */
if ictst .ne 1; goto loglik; endif;
@----------------------------------------
Johansen's Theorem II:
 Likelihood ratio test of H1: beta=Ah*phi
----------------------------------------@
hrk= rk*ah;
{beta2,lmbd2} = canon(hrk,r0,r,p);

@ test statistics @
tst2= nperg*sumc( ln( (1- lmbd2[1:r,1])./(1- lmbd1[1:r,1]) ));
df = r*(p-cols(ah))-ep; format /RDN 3,0;
pv = 100*(cdfchic(tst2,df));
" ";"LR test statistic of H1(beta=H*phi) when H ="; format /RDN 8,4; ah;
" ";"   Test statistic: ";;format /RDN 7,3; tst2;;
"       D.F.= ";;format /RDN 3,0;df;;
"       P-value(%)= ";; format /RDN 4,2; pv;format /RDN 11,5;
loglik:
@------------ Log-likelihoods ------------@

 @ form S00 @
S00=(r0'r0)/rows(r0);
@ Unconstrained log-likelihood @
ullik=-(nperg/2)*ln(det(s00))-(nperg/2)*sumc(ln(1-lmbd1));
" ";"------------ Log-likelihoods ------------";
"Unconstrained log-likelihood = ";;format 10,2;ullik;

@ Constrained log-likelihoods @
"Log-likelihood, when #of civ = 0, =";;-(nperg/2)*ln(det(s00));
cllik1=-(nperg/2)*ln(det(s00))-(nperg/2)*sumc(ln(1-lmbd1[1:r,1]));
"Log-likelihood, when # of civ is ";;format 1,0;r;;", = ";;
format 10,2; cllik1;
if ictst==1;
@ Subspace constrained Log-likelihood @
 cllik2=-(nperg/2)*ln(det(s00))-(nperg/2)*sumc(ln(1-lmbd2[1:r,1]));
 "Log-likelihood, when # of civ is ";;format /RDN 1,0;r;;format 10,2;
 " and H is given above, = ";;cllik2;
endif;
" ";
retp(beta);
endp;
