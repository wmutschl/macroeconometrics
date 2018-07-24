'Requires galiqje.wf1

pageselect gali_92

smpl 1955q1 1987q3

equation eq1.tsls  ygr ddrate dec1 dec2 ygr(-1 to -4) ddrate(-1 to -3) dec1(-1 to -3) dec2(-1 to -3) c @ c ygr(-1 to -4) drate(-4) ddrate(-1 to -3) ec1(-4) dec1(-1 to -3) dec2(-1 to -3) ec2(-4) 
eq1.makeresids eps1

equation eqvar.ls  ygr c ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
eqvar.makeresids res1

equation eq2.tsls  drate ygr diffec ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4) c @  c res1 eps1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
eq2.makeresids eps2

equation eq3.tsls  ec1 ygr drate ec2 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4) c @  c eps1 eps2 res1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
eq3.makeresids eps3

equation eq4.tsls ec2 ygr drate ec1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4) c @   c eps1 eps2 eps3 ygr(-1 to -4) drate(-1 to -4) ec1(-4) ec1(-1 to -4) ec2(-1 to -4)
 var galiqje.ls 1 4 ygr drate ec1 ec2  @ c 

scalar ca1=eq1.@coefs(1)
scalar ca2=eq1.@coefs(2)
scalar ca3=eq1.@coefs(3)

scalar ca4=eq2.@coefs(1)
scalar ca5=eq2.@coefs(2)

scalar ca6=eq3.@coefs(1)
scalar ca7=eq3.@coefs(2)
scalar ca8=eq3.@coefs(3)

scalar ca9=eq4.@coefs(1)
scalar ca10=eq4.@coefs(2)
scalar ca11=eq4.@coefs(3)

galiqje.cleartext(svar)

galiqje.append(svar) @e1=ca1*@e2+ca2*@e3+ca3*@e4+c(1)*@u1
galiqje.append(svar) @e2=ca4*@e1+ca5*@e3-ca5*@e4+c(2)*@u2
galiqje.append(svar) @e3=ca6*@e1+ca7*@e2+ca8*@e4+c(3)*@u3
galiqje.append(svar) @e4=ca9*@e1+ca10*@e2+ca11*@e3+c(4)*@u4

' f0=u means that one draws start values from a unifrom density , n=normal, 
galiqje.svar(rtype=text, f0=n,c=1e-08)
galiqje.results

'compute accumulated impulses 
galiqje.impulse(36, a, imp=struct, se=a)


