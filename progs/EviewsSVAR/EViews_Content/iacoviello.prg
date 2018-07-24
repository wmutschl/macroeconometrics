'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''IACOVIELLO (AER 2005) House prices, borrowing constraints, and monetary policy
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'wfclose "U:\My Documents\CURSOS\MFORE\OPR\Iacoviello\iacoviello.WF1"
'wfopen  "U:\My Documents\CURSOS\MFORE\OPR\Iacoviello\iacoviello.WF1"
'wfopen  "var_2014_bk_opr.WF1"


''It selects this page of the EVIEWS workfile to execute the instructions.

pageselect   Iacoviello
smpl 1970q1 2003q4

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''''''Here we transform the variables
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
genr fedbegin400 = ffr/400

genr LGDPDEF = LOG(GDPDEF)  
genr dLGDPDEF = d(lgdpdef)

genr LGDP = LOG(GDP) 
genr dLGDP = d(lgdp)

genr LHPI =  LOG(HPI) - LOG(GDPDEF)    
genr dlhpi = d(lhpi) 

genr lcrb= LOG( CRBSPOT)
genr dlcrb= d(lcrb) 

genr LPCE = LOG(PCE)
 

genr  Y    = 100*LGDP_F
genr  P    = 100*DLGDPDEF
genr  Q   = 100*LHPI_F
genr  R   = 100*FEDBEGIN400
genr CE = 100*LPCE_F

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''FIGURE 1 (IACOVIELLO) VAR EVIDENCE
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

smpl 1974q1 2003q2
!hrz=20
var varhousep1.ls  1  2  r p q y   @ c @trend LCRB(-1) D794AFTER
freeze(figa) varhousep1.impulse(!hrz, smat=ashk, se=mc,rep=400, matbys=a1)  r p q y @ r p q y
stop

'''''''''''Declares 16 vectors and fills them with the accumulated responses

for !i = 1 to 16

vector(!hrz) va{!i}    =  @columnextract(a1,{!i})
vector(!hrz) vaup{!i}    =  @columnextract(a1,{!i}) + 1.64485*(@columnextract(a1_se,{!i}))
vector(!hrz) valow{!i}   =  @columnextract(a1,{!i}) - 1.64485*@columnextract(a1_se,{!i})

next

''''''''''''''''builds matrices with 3 columns: the response and the bands
for !i= 1 to 16

     matrix(!hrz,3) g{!i} 	
	colplace( g{!i},  va{!i}, 1)
	colplace( g{!i},  vaup{!i}, 2)
	colplace( g{!i},  valow{!i}, 3)
'''''''''''''builds the individual charts
      freeze( fig!i) g{!i}.line
      fig{!i}.setelem(1)  lcolor(blue) lwidth(5) 
      fig{!i}.setelem(2)  lcolor(green) lwidth(3) 
      fig{!i}.setelem(3)  lcolor(green) lwidth(3) 
      fig{!i}.axis zeroline font(30)
      fig{!i}.options size(5,4) indenth(10.1)
      fig{!i}.legend -display   

next

 fig1.axis(left)  range(-0.2, 0.3)
 fig2.axis(left)  range(-0.2, 0.3)
 fig3.axis(left)  range(-1, 2)
 fig4.axis(left)  range(-1, 2)
 fig5.axis(left)  range(-0.2, 0.3)
 fig6.axis(left)  range(-0.2, 0.3)
 fig7.axis(left)  range(-1, 2)
 fig8.axis(left)  range(-1, 2)
 fig9.axis(left)  range(-0.2, 0.3)
 fig10.axis(left)  range(-0.2, 0.3)
 fig11.axis(left)  range(-1, 2)
 fig12.axis(left)  range(-1, 2)
 fig13.axis(left)  range(-0.2, 0.3)
 fig14.axis(left)  range(-0.2, 0.3)
 fig15.axis(left)  range(-1, 2)
 fig16.axis(left)  range(-1, 2)

 fig1.addtext(0.6, -0.9,,font(45,b))    response of R
 fig2.addtext(0.6, -0.9,,font(45,b))    response of PI
 fig3.addtext(0.6, -0.9, font(45,b))    response of q
 fig4.addtext(0.5, -0.9, font(45,b))    response of Y

fig1.addtext(l,b, font(45,b) )    "shock to R" 
fig5.addtext(l, font(45,b) )       "shock to PI"
fig9.addtext(l, font(45,b) )       "shock to q"
fig13.addtext(l, font(45,b) )     "shock to Y"

 fig13.addtext(15,b,font(45,b))  quarters
 fig14.addtext(15,b,font(45,b))  quarters
 fig15.addtext(15,b,font(45,b))  quarters
 fig16.addtext(15,b,font(45,b))  quarters

freeze(figR) fig1 fig2 fig3 fig4 fig5 fig6  fig7 fig8 fig9 fig10 fig11 fig12 fig13 fig14 fig15 fig16
figR.align(4,2,1.0)
figR.options  indenth(90.1,5)
figR.addtext(8.0,20.9, font(40))  Figure I.  VAR Evidence, UNITED STATES
figR.addtext(-0.0,21.7, font(35)) Notes: VAR estimated from 1974Q1 to 2003Q2. The dashed lines indicate 90-percent confidence bands. The Choleski
figR.addtext(-0.0,22.3, font(35)) ordering of the impulse responses is R, PI, q, Y. Coordinate: percent deviation from the baseline. 
show figR

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
'''''''''''''FIGURE 3 (IACOVIELLO) RESPONSE OF AGGREGATE CONSUMPTION TO HOUSING PRICE SHOCK
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
var varhousep2.ls  1  2  r q ce y p   @ c @trend LCRB(-1) D794AFTER

 
freeze(figd) varhousep2.ls  1  2  r q ce y p   @ c @trend LCRB(-1) D794AFTER


 
!hrz=10
freeze(figb) varhousep2.impulse(!hrz,  se=mc,rep=400, matbys=a2)  r q ce y p  

''''''''''''''''''''''''''''''''''''''
'''''''''''A2 and A2_SE are  the matrices obtained with the estimated VAR 
'''''''''''Declares 16 vectors and fills them with the accumulated responses
'''''''''''''''''''''''''''''''''''''
for !i = 1 to  25

vector(!hrz) vaa{!i}           =  @columnextract(a2,{!i})
vector(!hrz) vaupa{!i}      =  @columnextract(a2,{!i}) + 1.64485*(@columnextract(a2_se,{!i}))
vector(!hrz) valowa{!i}     =  @columnextract(a2,{!i}) - 1.64485*@columnextract(a2_se,{!i})

next


for !i= 1 to 25
     matrix(!hrz,3) ga{!i} 	
	colplace( ga{!i},  vaa{!i}, 1)
	colplace( ga{!i},  vaupa{!i}, 2)
	colplace( ga{!i},  valowa{!i}, 3)
next

'''''''''''''''Here we normalize the responses by the variance of the shock
'''''''''''''''Equivalente to normalizing Matrix B.

      matrix(!hrz,3) gaa7 = ga7/ga7(1,1)
      matrix(!hrz,3) gaa8 = ga8/ga7(1,1)
 
''''''''''''''Here we build charts 7 and 8
for !i = 7 to 8
      freeze( figaa!i) gaa{!i}.line
      figaa{!i}.setelem(1)  lcolor(blue) lwidth(3) 
      figaa{!i}.setelem(2)  lcolor(green) lwidth(3) 
      figaa{!i}.setelem(3)  lcolor(green) lwidth(3) 
      figaa{!i}.axis zeroline font(30)
      figaa{!i}.options size(15,5) 
      figaa{!i}.legend -display   
 next
 
    figaa7.axis(left)  range(-1, 2)
    figaa8.axis(left)  range(-0.5, 0.8)
    figaa7.addtext(4.7, -0.9,font(45,b))    HOUSE PRICES
    figaa8.addtext(4.7, -0.9, font(45,b))    CONSUMPTION

''''''''''This time we will build the figure putting together only charts 7 and 8

   freeze(figaR3) figaa7 figaa8


   figaR3.addtext(-.5,0.7,5,7 ,b,font(40))   FIGURE 3. Response of aggregate consumption to a housing price
   figaR3.addtext(3.2,13.2 ,font(40))  shock: various values of m and m"
   figaR3.addtext(0.7, 6,5 font(95))   
   show figaR3


