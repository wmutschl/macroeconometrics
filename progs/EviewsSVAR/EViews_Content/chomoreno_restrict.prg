'Requires chomoreno.wf1
smpl 1981q1 2001q1

equation eqdgp.ls gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c
eqdgp.makeresids res1

equation eq2.ls infl gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c
eq2.makeresids eps2

equation eq3.tsls ff gap infl gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c @ res1 eps2 gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)
eq3.makeresids eps3

equation eq1.tsls gap infl ff gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c @ eps2 eps3 gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)

scalar ca1=eq1.@coefs(1)
scalar ca3=eq3.@coefs(1)
scalar ca4=eq3.@coefs(2)

var chomorrest.ls 1 2 gap infl ff

chomorrest.cleartext(svar)
chomorrest.append(svar) @e1=ca1*@e2 + c(1)*@u1
chomorrest.append(svar) @e2=c(2)*@u2
chomorrest.append(svar) @e3=ca3*@e1+ca4*@e2+c(3)*@u3
' f0=u means that one draws start values from a unifrom density , n=normal, 
chomorrest.svar(rtype=text, f0=n)
'compute normal impulses
chomorrest.impulse(10,imp=struct, se=a)


