'Requires bqdata.wf1

smpl 1950q2 1987q4

equation eq1.tsls dya du dya(-1 to -8) du(-1 to -7)  c @ dya(-1 to -8) du(-1 to -7) u(-8)
eq1.makeresids eps1

equation eq2.tsls u dya dya(-1 to -8) u(-1 to -8) c @ dya(-1 to -8) u(-1 to -8) eps1
eq2.results

scalar ca=eq1.@coefs(1)
scalar cb=eq2.@coefs(1)

var bq.ls 1 8  dya u @ c
bq.cleartext(svar)
bq.append(svar) @e1=ca*@e2+c(1)*@u1
bq.append(svar) @e2=cb*@e1+c(2)*@u2

' f0=u means that one draws start values from a uniform density , n=normal, 
bq.svar(rtype=text, f0=u)
'compute normal impulses
bq.impulse(36,imp=struct, se=a)

