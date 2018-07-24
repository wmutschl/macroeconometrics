'Requires gdp_m2.wf1
equation eq1.tsls s_dgdp d(s_dm2)  s_dgdp(-1)  c @ s_dgdp(-1) s_dm2(-1)
eq1.makeresids eps1

equation eq2.tsls s_dm2 s_dgdp s_dgdp(-1) s_dm2(-1) c @ s_dgdp(-1) s_dm2(-1) eps1

scalar ca=eq1.@coefs(1)
scalar cb=eq2.@coefs(1)

var gdpm.ls 1 1  s_dgdp s_dm2

gdpm.cleartext(svar)
gdpm.append(svar) @e1=ca*@e2+c(1)*@u1
gdpm.append(svar) @e2=cb*@e1+c(2)*@u2

gdpm.svar(rtype=text,f0=u)
'compute normal impulses
gdpm.impulse(36,imp=struct, se=a)
'compute accumulated impulses 
gdpm.impulse(36, a, imp=struct, se=a)
'compute accumulated impulses 
gdpm.impulse(36, a, imp=struct, se=a) @ 1 2
'

