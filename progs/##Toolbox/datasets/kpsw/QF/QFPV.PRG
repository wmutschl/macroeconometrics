loadp eig, datestr, timestr;

cf=files(ifname,0); if cf[1,1] $/= 0;
if trunc((maxl-minl)/incl)*incl+minl .eq maxl;

if (vals(org) .eq  79) or (vals(org) .eq 111);
 load d0[altim,ns] = ^ifname; d2=d0;
elseif (vals(org) .eq 86) or (vals(org) .eq 118);
 load d0[ns,altim] = ^ifname; d2=d0'; endif;
d1=d2[bg:ed,.]; clear d0,d2;
d=d1;  clear d1; @ d contains the notim,ns data vector @

@ Normalizing data @
notim=ed-bg+1;

if (optn .eq 1); sfname = "\\kpsw\\qf\\sort" $+ ftocv(k,1,0) $+ "t";
@ regress on one and t @
xt=ones(notim,1)~seqa(1,1,notim); dbeta=solpd(xt'd,xt'xt);
dn=d-xt*dbeta;

elseif (optn .eq 2); sfname = "\\kpsw\\qf\\sort" $+ ftocv(k,1,0) $+ "c";
@ regress on one @
xt=ones(notim,1); dbeta=solpd(xt'd,xt'xt);
dn=d-xt*dbeta;

elseif (optn .eq 3); sfname = "\\kpsw\\qf\\sort" $+ ftocv(k,1,0);
@ do nothing @
dn=d;
endif;  clear d, xt;

cf=files(sfname,0); if cf[1,1] $/= 0;
if (notim*maxl*k)>2e+4; lflag=1; else; lflag=0; endif;

@ form estimate of phi using all of the data to get CI vectors @
ya=dn;

@ Estimate phi @
zd=ya[1:notim-1,.]; ye=ya[2:notim,.]; phi=solpd(zd'ye,zd'zd);

clear evr,evi;
ev=zeros(ns,ns); ee=eig(phi,"ev"); evalm=sqrt((ee[.,1].^2)+(ee[.,2].^2));

@ Check whether correct EIG proc is loaded @
if rows(evr) == rows(ee);

@ Form Potential CI Vectors @
nc=1; do while nc <= ns;
if ee[nc,2] /= 0;
ev[.,nc]=evr[.,nc]+evi[.,nc]; ev[.,nc+1]=evr[.,nc]-evi[.,nc];
nc=nc+2;
else;
ev[.,nc]=evr[.,nc];
nc=nc+1;
endif;
endo;

ee=evalm~ev'; es=sortc(ee,1); ee=es';
ev[1:ns,1:ns]=ee[2:ns+1,1:ns];

@ form AR - Principle components -- Z@
z=dn*ev;     clear dn, ev, evi, evr;

nl=0;
ya=z[1+nl:notim,.];

@ Estimate phi using last k filtered values of z @
x=ya[1:notim-1-nl,ns-k+1:ns]; y=ya[2:notim-nl,ns-k+1:ns];
phi=solpd(x'y,x'x);

ee=eig(phi,0); es=sortc(ee,1);
qf0=notim*(es[.,1]-1);

@ Difference the last k elements of Z using this estimates PHI @
zd=z[2:notim,ns-k+1:ns]-z[1:notim-1,ns-k+1:ns]*phi;

@ form first estimate of phi using nl = 0 @

@ Form autoregressions of e @
nl=minl; do while nl <= maxl;
nl1=1; do while nl1 <= nl;
xe=zd[1+nl-nl1:notim-1-nl1,.];
oxe = "xe" $+ ftocv(nl1,2,0);

if lflag; save ^oxe=xe;
 else; vp=varput(xe,oxe);
endif; clear xe;

nl1=nl1+1; endo;
   dumyr=zeros(1,k);
 jj=1;
 do while jj <= nl; xexe=zeros(nl*k,nl*k);
    ii=1;
    do while ii <= jj;
    if lflag;
      nxeii = "xe" $+ ftocv(ii,2,0);
      nxejj = "xe" $+ ftocv(jj,2,0);
      load xeii=^nxeii;
      load xejj=^nxejj;
    else;
      xeii = varget("xe" $+ ftocv(ii,2,0));
      xejj = varget("xe" $+ ftocv(jj,2,0));
    endif;
      xeji = xejj'xeii;
      oxeji="xe" $+ ftocv(jj,2,0) $+ ftocv(ii,2,0);
      vp=varput(xeji,oxeji);
        if jj .ne ii;
          oxeij = "xe" $+ ftocv(ii,2,0) $+ ftocv(jj,2,0);
          xeij = xeji';
          vp=varput(xeij,oxeij);
        endif;
    ii = ii+1;
    endo;
  jr=1; do while jr<=jj;
  ir=1; do while ir<=jj;
    xexe[(jr-1)*k+1:jr*k,(ir-1)*k+1:ir*k] =
            varget("xe"$+ftocv(jr,2,0)$+ftocv(ir,2,0));
  ir=ir+1; endo;
  jr=jr+1; endo;
 ye=zd[1+nl:notim-1,.];
 xye = xeii'ye;
 dumyr=dumyr | xye;
 jj = jj+1;
 endo;         clear xye, xeii, xejj, xeij, ye;
xye = dumyr[2:nl*k+1,.];
rho=solpd(xye,xexe); clear xexe, dumyr;

@ filter z's @
ya=z[1+nl:notim,ns-k+1:ns];
nl1=1; do while nl1 <= nl;
yb=ya;
ns1=(nl1-1)*k;
yfit=z[1+nl-nl1:notim-nl1,ns-k+1:ns]*rho[ns1+1:ns1+k,1:k];
ya=yb-yfit;
nl1=nl1+1;endo;   clear yb, yfit;

@ Estimate phi using last k filtered values of z@

x=ya[1:notim-1-nl,.]; y=ya[2:notim-nl,.];
phi=solpd(x'y,x'x);

ee=eig(phi,0); es=sortc(ee,1);
qf=notim*(es[.,1]-1);
qf0=qf0~qf;
nl=nl+incl;
endo;         clear z, zd, ya, x, y, xye;

@ Load test stat table. @
load srt[399,k+1] = ^sfname;

pv=zeros(k,6);
i=1; do while i<=k;
        pv[i,1] = counts(srt[.,i+1] , qf0[i,1])./4.0;
   j=2; do while j<= cols(qf0);
        pv[i,j] = counts(srt[.,i+1] , qf0[i,j])./4.0;
   j=j+1;
   endo;
i=i+1;
endo;    clear srt;

@ Procedure to print qf and p-value. @
proc(0) = qfpv(i); local jp,jpx;
" ";
"      Number of lags                QF              P-value(%)";
"      ==============            ==========          ==========";
"            ";;"0";;"                    ";;
format /rd 8,4;qf0[i,1];;
"           ";;format /rd 6,2;pv[i,1];
jp=2; do while jp<=2+(maxl-minl)/incl;jpx=minl+(jp-2)*incl;
"           ";;format /rd 2,0; jpx;;"                   ";;
format /rd 8,4;qf0[i,jp];;
"           ";;format /rd 6,2;pv[i,jp];
jp=jp+1;
endo;
" "; " "; " "; retp; endp;

output file = ^ofname reset; @ screen off;@
@ output on; @
format /rz 5,0;
"DATE : ";;datestr; "TIME : ";;timestr;
" "; " "; " "; " ";
"                                 Data file = ";; upper(ifname);
"      Total number of time periods in file =";; altim;
"                First period used for test =";; bg;
"                 Last period used for test =";; ed;
"                     Number of data series =";; ns;
"                             Preprocessing = ";;
if (optn .eq 1); "Linear detrending"; elseif (optn .eq 2); "Mean removed";
elseif (optn .eq 3); "None"; endif;
"                 H0 : Number of unit roots =";; k;
" "; " "; " ";
h=1; do while h<=k;
"      H";;format 1,0;h;;": Number of unit roots = ";; kh=k-h; kh;
 call qfpv(h); h=h+1; endo;
@output off; screen on;@

if lflag; dos del xe*.fmt; endif;

else;
" ERROR: THE PROGRAM CANNOT ACCESS EIGEN VECTORS. ";
"        PLEASE CHECK WHETHER EIG PROCEDURE CAN GENERATE EIGEN VECTORS, ";
"        WHICH SHOULD BE GLOBAL SYMBOLS."; endif;
else; "ERROR: THE FILE " upper(sfname) " CANNOT BE FOUND."; endif;
else; "ERROR: THE LAG INCREMENT IS NOT CORRECT."; endif;
else; "ERROR: THE FILE " upper(ifname) " CANNOT BE FOUND."; endif;

