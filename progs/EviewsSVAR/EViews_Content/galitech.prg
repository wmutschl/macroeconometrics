'Requires galidusa.wf1
'Replicate Gali (1999) using the IV approach.

pageselect orig_gali

smpl @all

genr lprodh = log(gdpq)-log(lpmhu)
genr dp = 100*dlog(p) 'Inflation
genr dinf = d(dp)
genr ddinf = d(dinf)
genr dprodh = 100*d(lprodh)
genr dhours = 100*dlog(lpmhu)
genr ddhours = d(dhours)
genr dm2 = 100*dlog(m2)
genr ec1 = (rm3/4)-dp
genr ec2 = dm2 - dp
genr dec1 = d(ec1)
genr dec2 = d(ec2)

smpl 1959Q1 1994Q4

equation eq1.tsls  dprodh  ddhours ddinf dec1 dec2 ddhours(-1 to -3) ddinf(-1 to -3)  dec1(-1 to -3)  dec2(-1 to -3) dprodh(-1 to -4) c @ c dprodh(-1 to -4) dhours(-1 to -4) dinf(-1 to -4) ec1(-4) dec1(-1 to -3) ec2(-4) dec2(-1 to -3)
show eq1.results

var galitech.ls 1 4 dprodh dhours dinf ec1 ec2  @ c 

scalar ca1=eq1.@coefs(1)
scalar ca2=eq1.@coefs(2)
scalar ca3=eq1.@coefs(3)
scalar ca4=eq1.@coefs(4)

galitech.cleartext(svar)

galitech.append(svar) @e1=ca1*@e2+ca2*@e3 + ca3*@e4 +ca4*@e5+c(1)*@u1         
galitech.append(svar) @e2=c(2)*@u1+c(3)*@u2
galitech.append(svar) @e3= c(4)*@u1+c(5)*@u2+c(6)*@u3
galitech.append(svar) @e4= c(7)*@u1+c(8)*@u2+c(9)*@u3+c(10)*@u4
galitech.append(svar) @e5= c(11)*@u1+c(12)*@u2+c(13)*@u3+c(14)*@u4+c(15)*@u5

galitech.svar(rtype=text, f0=n)
galitech.results
'compute accumulated impulses 
galitech.impulse(12, a, imp=struct, se=a)

