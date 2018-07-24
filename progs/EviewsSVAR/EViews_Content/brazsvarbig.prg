

scalar vv

 scalar nv
nv=7


'computing impulse responses

matrix (nv,nv) a0
matrix (nv,nv) a1



matrix (nv,nv) b

' define A0 in terms of the C matrix
a0(1,1)=1
a0(2,1)=-C(55)
a0(2,2)=1
a0(3,1)=-C(9)
a0(3,2)=-C(10)
a0(3,3)=1
a0(4,1)=-C(19)
a0(4,2)=-C(20)
a0(4,3)=-C(57)
a0(4,4)=1
a0(5,1)=-C(29)
a0(5,2)=-C(30)
a0(5,3)=-C(58)
a0(5,4)=-C(59)
a0(5,5)=1
A0(6,1)=-C(39)
A0(6,2)=-C(40)
A0(6,3)=-C(60)
A0(6,4)=-C(61)
A0(6,5)=-C(62)
A0(6,6)=1

A0(7,1)=-C(49)
A0(7,2)=-C(50)
A0(7,3)=-C(63)
A0(7,4)=-C(64)
A0(7,5)=-C(65)
A0(7,6)=-C(66)
A0(7,7)=1

a1(1,1)=C(53)
a1(1,2)=C(54)
a1(2,1)=C(56)
a1(2,2)=C(57)

'construct A1 and A2 from the C matrix
vv=0
 '3rd equation
for !j=3 to nv
a1(3,!j) = C(!j+vv-2)
next 

A1(3,1)=C(7)
A1(3,2)=C(8)

'4th equation
vv=10
for !j=3 to nv
a1(4,!j) = C(!j+vv-2)
next
A1(4,1)=C(17)
A1(4,2)=C(18)

 '5th equation
vv=20
for !j=3 to nv
a1(5,!j) = C(!j+vv-2)
next 
a1(5,1)=C(27)
a1(5,2)=C(28)

'6th equation
vv=30
for !j=3 to nv
a1(6,!j) = C(!j+vv-2)
next 
a1(6,1)=C(37)
a1(6,2)=C(38)
'7th equation
vv=40
for !j=3 to nv
a1(7,!j) = C(!j+vv-2)
next 
a1(7,1)=C(47)
a1(7,2)=C(48)

'if want unit shocks put diag(b) to unity. If std dev shocks put to that 
b(1,1)=.409
b(2,2)=1
b(3,3)=1
b(4,4)=1
b(5,5)=1
b(6,6)=1
b(7,7)=1





matrix(nv,nv)b1
'compute the VAR coefficients from the SVAR coefficients

b1=@inverse(a0)*a1


'set impulse response horizon

scalar horz
horz=100

matrix (nv,nv) cc0


'form initial impulse responses


cc0=(@inverse(a0))*b


matrix(horz,nv) imp1
matrix(horz,nv) imp2
matrix(horz,nv) imp3
matrix(horz,nv) imp4
matrix(horz,nv) imp5
matrix(horz,nv) imp6
matrix(horz,nv) imp7




'this stuff must be able to go into a subroutine
' it just extracts the impulse responses for each shock from the matrices

vector v1=@columnextract(cc0,1)
vector v2=@columnextract(cc0,2)
vector v3=@columnextract(cc0,3)
vector v4=@columnextract(cc0,4)
vector v5=@columnextract(cc0,5)
vector v6=@columnextract(cc0,6)
vector v7=@columnextract(cc0,7)



for !j=1 to nv
imp1(1,!j) = v1(!j)
imp2(1,!j) = v2(!j)
imp3(1,!j) = v3(!j)
imp4(1,!j) = v4(!j)
imp5(1,!j) = v5(!j)
imp6(1,!j) = v6(!j)
imp7(1,!j) = v7(!j)
next 
'now compute 1 and > period ahead impulse responses. These follow C1=b1*C0
matrix(nv,nv) cc1

for !kk=2 to horz

cc1=b1*cc0

vector v1=@columnextract(cc1,1)
vector v2=@columnextract(cc1,2)
vector v3=@columnextract(cc1,3)
vector v4=@columnextract(cc1,4)
vector v5=@columnextract(cc1,5)
vector v6=@columnextract(cc1,6)
vector v7=@columnextract(cc1,7)

for !j=1 to nv
imp1(!kk,!j) = v1(!j)
imp2(!kk,!j) = v2(!j)
imp3(!kk,!j) = v3(!j)
imp4(!kk,!j) = v4(!j)
imp5(!kk,!j) = v5(!j)
imp6(!kk,!j) = v6(!j)
imp7(!kk,!j) = v7(!j)
next 
cc0=cc1
next



'imp1 shows effects of first shock on the nv variables  The firth shock deson't exist due to the identity
'rows are periods and columns are the variables
show imp1

