'Requires galiqje.wf1
pageselect gali_alt
smpl 1955Q1 1987Q3

equation alt_eq1.tsls ygr ddrate dec1 dec2 ygr(-1 to -4) ddrate(-1 to -3) dec1(-1 to -3) dec2(-1 to -3) c @ c ygr(-1 to -4) drate(-4) ddrate(-1 to -3) ec1(-4) dec1(-1 to -3) dec2(-1 to -3) ec2(-4) 
alt_eq1.makeresids eps1

equation alt_eqvar.ls ygr c ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
alt_eqvar.makeresids res1

equation alt_eq2.tsls drate ygr dec1 dec2  ygr(-1 to -4) drate(-1 to -4) dec1(-1 to -3) dec2(-1 to -3) c @  c  eps1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
alt_eq2.makeresids eps2lr

equation alt_eq3.tsls ec1 ygr drate ec2 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4) c @  c eps1 eps2lr res1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4)
alt_eq3.makeresids eps3

equation alt_eq4.tsls ec2 ygr drate ec1 ygr(-1 to -4) drate(-1 to -4) ec1(-1 to -4) ec2(-1 to -4) c @   c eps1 eps2lr eps3 ygr(-1 to -4) drate(-1 to -4) ec1(-4) ec1(-1 to -4) ec2(-1 to -4)

var galiqjealt.ls 1 4 ygr drate ec1 ec2  @ c

scalar alt_ca1=alt_eq1.@coefs(1)
scalar alt_ca2=alt_eq1.@coefs(2)
scalar alt_ca3=alt_eq1.@coefs(3)

scalar alt_ca4=alt_eq2.@coefs(1)
scalar alt_ca5=alt_eq2.@coefs(2)
scalar alt_ca6=alt_eq2.@coefs(3)

scalar alt_ca7=alt_eq3.@coefs(1)
scalar alt_ca8=alt_eq3.@coefs(2)
scalar alt_ca9=alt_eq3.@coefs(3)

scalar alt_ca10=alt_eq4.@coefs(1)
scalar alt_ca11=alt_eq4.@coefs(2)
scalar alt_ca12=alt_eq4.@coefs(3)

galiqjealt.cleartext(svar)

galiqjealt.append(svar) @e1=alt_ca1*@e2+alt_ca2*@e3+alt_ca3*@e4+c(1)*@u1
galiqjealt.append(svar) @e2=alt_ca4*@e1+alt_ca5*@e3+alt_ca6*@e4+c(2)*@u2
galiqjealt.append(svar) @e3=alt_ca7*@e1+alt_ca8*@e2+alt_ca9*@e4+c(3)*@u3
galiqjealt.append(svar) @e4=alt_ca10*@e1+alt_ca11*@e2+alt_ca12*@e3+c(4)*@u4

' f0=u means that one draws start values from a unifrom density , n=normal, 
galiqjealt.svar(rtype=text, f0=n)
galiqjealt.results
'compute accumulated impulses 
galiqjealt.impulse(36, a, imp=struct, se=a)


