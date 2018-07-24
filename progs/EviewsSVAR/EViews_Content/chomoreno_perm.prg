'
'Requires chomoreno.wf1
smpl 1981q1 2001q1
equation eq1.tsls dgap dff dinfl dff(-1) dgap(-1 to -2) dinfl( -1) c @ c dgap(-1 to -2) infl(-2) dinfl(-1) ff(-2) dff(-1) 

eq1.makeresids eps1
equation eqdgp.ls  dgap c dgap(-1 to -2) infl(-1 to -2) ff(-1 to -2)
eqdgp.makeresids res1

equation eq2.tsls ff dgap infl dgap(-1 to -2) infl(-1 to -2) ff(-1 to -2) c @ c eps1 res1 dgap(-1 to -2) infl( -1 to -2) ff(-1 to -2)
eq2.makeresids eps2
 

equation eq3.tsls infl dgap ff dgap(-1 to -2) infl(-1 to -2) ff(-1 to -2) c @  c eps1 eps2 dgap(-1 to -2) infl( -1 to -2) ff(-1 to -2)
eq3.makeresids eps3

var chomorperm.ls 1 2 dgap ff infl

scalar ca1=eq1.@coefs(1)  'ff 
scalar ca2=eq1.@coefs(2)  'inf

scalar ca3=eq2.@coefs(1) 'gap
scalar ca4=eq2.@coefs(2) 'inf

scalar ca5=eq3.@coefs(1) 'gap
scalar ca6=eq3.@coefs(2) 'ff

chomorperm.cleartext(svar)

chomorperm.append(svar) @e1=ca1*@e2+ca2*@e3+c(1)*@u1
chomorperm.append(svar) @e2=ca3*@e1+ca4*@e3-+c(2)*@u2
chomorperm.append(svar) @e3=ca5*@e1+ca6*@e2+c(3)*@u3


' f0=u means that one draws start values from a unifrom density , n=normal, 

chomorperm.svar(rtype=text, f0=n,c=1e-10)

'compute normal impulses
chomorperm.impulse(36,imp=struct, se=a)


