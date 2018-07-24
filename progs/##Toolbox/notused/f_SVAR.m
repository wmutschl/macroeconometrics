function f=f_SVAR(B0inv,SIGMAUHAT,A1inv,Rshort,Rlong,selB0inv,selTheta)
Theta = A1inv*B0inv;

f=[vech(SIGMAUHAT-B0inv*B0inv');
   B0inv(selB0inv) - Rshort(selB0inv);
   Theta(selTheta) - Rlong(selTheta);
   ];
    