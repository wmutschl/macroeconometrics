
'Computing Impulse Responses for Restricted Brazil VAR
' nvar= no endogenous variables
' nex = number exogenous variables ( not constant, dummies or deterministic trends
'lag = order VAR
'horz= horizon for impulse responses


 scalar nvar
scalar nex
scalar nx
scalar lag
scalar horz
lag=1
horz=40

nvar=5

nex=2
nx=nvar+nex



' need to compute B1 in z(t)=B1*z(t-1)+F*zi(t)

'zi has the shocks first, exog vars last
'there are nx of them


matrix (nvar,nvar) b1






'construct B1 and F from the C matrix

' need to set vv to be 1 less that first C(..) in the equation being described

scalar vv
vv=0
 'first equation
for !j=1 to nvar
b1(1,!j) = C(!j+vv)
next 
'2nd equation
vv=7
for !j=1 to nvar
b1(2,!j) = C(!j+vv)
next
 '3rd equation
vv=14
for !j=1 to nvar
b1(3,!j) = C(!j+vv)
next 

'4th equation
vv=21
for !j=1 to nvar
b1(4,!j) = C(!j+vv)
next 

'5th equation
vv=28

for !j=1 to nvar
b1(5,!j) = C(!j+vv)
next 

matrix(nvar,nx) f

' this will produce 1 unit shocks. If want std dev need to multiply these elements by std deviations
' std dev of errors are 1.084, .492, 3.552, 1.349, 8.441
'std dev exog vars is 1.28 (ystar) and 1.676 (rus)

f(1,1)=1
f(2,2)=1
f(3,3)=1
f(4,4)=1
f(5,5)=1
f(1,6)=C(7)
f(2,6)=C(14)
f(3,6)=C(21)
f(4,6)=C(28)
f(5,6)=C(35)
f(5,7)=C(36)

'set horizon for impulse responses

'compute contempraneous responses

matrix (nvar,nex) cc0


'form initial impulse responses

cc0=f

'there are 5 shocks and 2 exog variable "shocks"
'so impj is response of the 5 variables to the seven shocks

matrix(horz,nvar) imp1
matrix(horz,nvar) imp2
matrix(horz,nvar) imp3
matrix(horz,nvar) imp4
matrix(horz,nvar) imp5
matrix(horz,nvar) imp6
matrix(horz,nvar) imp7



'this stuff must be able to go into a subroutine
' it just extracts the impulse responses for each shock from the matrices

vector v1=@columnextract(cc0,1)
vector v2=@columnextract(cc0,2)
vector v3=@columnextract(cc0,3)
vector v4=@columnextract(cc0,4)
vector v5=@columnextract(cc0,5)
vector v6=@columnextract(cc0,6)
vector v7=@columnextract(cc0,7)


for !j=1 to nvar
imp1(1,!j) = v1(!j)
imp2(1,!j) = v2(!j)
imp3(1,!j) = v3(!j)
imp4(1,!j) = v4(!j)
imp5(1,!j) = v5(!j)
imp6(1,!j) = v6(!j)
imp7(1,!j) = v7(!j)
next 

'now compute 1 period ahead and following impulse responses. These follow C1=b1*C0



matrix(nvar,nvar) cc1


for !kk=lag to horz
cc1=b1*cc0
vector v1=@columnextract(cc1,1)
vector v2=@columnextract(cc1,2)
vector v3=@columnextract(cc1,3)
vector v4=@columnextract(cc1,4)
vector v5=@columnextract(cc1,5)
vector v6=@columnextract(cc1,6)
vector v7=@columnextract(cc1,7)

for !j=1 to nvar
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


'imp1 shows effects of first shock on the 5 variables  
'rows are periods and columns are the variables
show imp1
