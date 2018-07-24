'Requires gdp_m2.wf1

gdpmiv.tsls
scalar ca1=gdpmiv.@coefs(2)
var gdpm.ls 1 1  s_dgdp s_dm2
gdpm.results
gdpm.cleartext(svar)
gdpm.append(svar) @e1=ca1*@e2+c(1)*@u1
gdpm.append(svar) @e2=c(3)*@e1+c(2)*@u2
gdpm.svar(rtype=text,f0=u)
'compute accumulated impulses 
gdpm.impulse(36, a, imp=struct, se=a)


