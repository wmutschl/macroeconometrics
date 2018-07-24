'Requires peersman.wf1
smpl 1980q1 2002q2
equation eq1.ls dpoil  dpoil(-1 to -3) dusgdp(-1 to -3) duscpi(-1 to -3) usint(-1 to -3) c  time
eq1.makeresids eps1

equation eq2.tsls dusgdp dusint dduscpi dpoil dpoil(-1 to -3) dusgdp(-1 to -3) dusint(-1 to -2) dduscpi(-1 to -2) c  time @  c eps1 dpoil(-1 to -3) dusgdp(-1 to -3)  duscpi(-1 to -3) usint(-1 to -3)  time
eq2.results
eq2.makeresids eps2
 
equation eqdus.ls  dusgdp c dpoil(-1 to -3) dusgdp(-1 to -3) usint(-1 to -3) duscpi(-1 to -3)  time
eqdus.makeresids res2

equation eq3.tsls usint  dpoil dusgdp duscpi dpoil(-1 to -3) dusgdp(-1 to -3) usint(-1 to -3) duscpi(-1 to -3)  c time @  eps1 eps2 res2 dpoil(-1 to -3) dusgdp(-1 to -3) usint(-1 to -3) duscpi(-1 to -3) c time
eq3.results
eq3.makeresids eps3

equation eq4.tsls duscpi dpoil dusgdp usint  dpoil(-1 to -3) dusgdp(-1 to -3) usint(-1 to -3) duscpi(-1 to -3)  c time @ eps1 eps2 eps3 dpoil(-1 to -3) dusgdp(-1 to -3) usint(-1 to -3) duscpi(-1 to -3)  c time
eq4.results

var peersman.ls 1 3 dpoil dusgdp usint duscpi   @ c time

scalar ca1=eq2.@coefs(3)
scalar ca2=eq2.@coefs(1)
scalar ca3=eq2.@coefs(2)

scalar ca4=eq3.@coefs(1)
scalar ca5=eq3.@coefs(2)
scalar ca6=eq3.@coefs(3)

scalar ca7=eq4.@coefs(1)
scalar ca8=eq4.@coefs(2)
scalar ca9=eq4.@coefs(3)

peersman.cleartext(svar)

peersman.append(svar) @e1=c(1)*@u1
peersman.append(svar) @e2=ca1*@e1+ca2*@e3+ca3*@e4+c(2)*@u2
peersman.append(svar) @e3=ca4*@e1+ca5*@e2+ca6*@e4+c(3)*@u3
peersman.append(svar) @e4=ca7*@e1+ca8*@e2+ca9*@e3+c(4)*@u4

' f0=u means that one draws start values from a unifrom density , n=normal, 
peersman.svar(rtype=text, f0=n,c=1e-10)
peersman.results
'compute normal impulses
'peersman.impulse(50,imp=struct, se=a)@ 1 2 3
peersman.impulse(30,imp=struct, se=a)
'compute accumulated impulses 
peersman.impulse(30, a, imp=struct, se=a)

