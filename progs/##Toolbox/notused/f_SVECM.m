function f=f_SVECM(B0inv,SIGMAUHAT,XI,Rshort,Rlong,selB0inv,selUpsilon)
UPSILON = XI*B0inv;

f=[vech(SIGMAUHAT-B0inv*B0inv');
   B0inv(selB0inv)     - Rshort(selB0inv);
   UPSILON(selUpsilon) - Rlong(selUpsilon);
   ];
    