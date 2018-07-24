smpl 1959:01 2001:7

scalar nvars = 120
scalar nslow = 70

!imp_periods = 48
!var_lags = 13  ''Number of lags is uncertain; they mention 13 in the paper

for %s ipp ipf ipc ipcd ipcn ipe ipi ipm ipmd ipmnd ipmfg ipd ipn ipmin iput ip _
ipxmca pmi pmp gmpyq gmyxpq lhel lhelx lhem lhnag lhur lhu680 lhu5 lhu14 _
lhu15 lhu26 lpnag lp lpgd lpmi lpcc lpem lped lpen lpsp lptu lpt lpfr lps lpgov _
lphrm lpmosa pmemp gmcq gmcdq gmcnq gmcsq gmcanq hsfr hsne hsmw _'
hssou hswst hsbr hmob pmnv pmno pmdel mocmq msondq fsncom fspcom _
fspin fspcap fsput fsdxp fspxe exrsw exrjan exruk exrcan fyff fygm3 fygm6 fygt1 _
fygt5 fygt10 fyaaac fybaac sfygm3 sfygm6 sfygt1 sfygt5 sfygt10 sfyaaac sfybaac _
fm1 fm2 fm3 fm2dq fmfba fmrra fmrnba fclnq fclbmc ccinrv pmcp pwfsa pwfcsa _ 
pwimsa pwcmsa psm99q punew pu83 pu84 pu85 puc pucd pus puxf puxhs _
puxm lehcc lehm hhsntn

 	genr sd_{%s} = ({%s} - @mean({%s}))/ @stdev({%s})
'	genr sd_{%s} = {%s} 'Uncomment this line to bypass standardization....

next

!dampener = 0.25/@stdev(fyff)  '25 basis points

'The first step is to extract principal components from the entire dataset (120 variables, including the federal funds rate)

group x_series_f3 sd_ipp sd_ipf sd_ipc sd_ipcd sd_ipcn sd_ipe sd_ipi sd_ipm sd_ipmd sd_ipmnd sd_ipmfg sd_ipd sd_ipn sd_ipmin sd_iput sd_ip _
sd_ipxmca sd_pmi sd_pmp sd_gmpyq sd_gmyxpq sd_lhel sd_lhelx sd_lhem sd_lhnag sd_lhur sd_lhu680 sd_lhu5 sd_lhu14 _
sd_lhu15 sd_lhu26 sd_lpnag sd_lp sd_lpgd sd_lpmi sd_lpcc sd_lpem sd_lped sd_lpen sd_lpsp sd_lptu sd_lpt sd_lpfr sd_lps sd_lpgov _
sd_lphrm sd_lpmosa sd_pmemp sd_gmcq sd_gmcdq sd_gmcnq sd_gmcsq sd_gmcanq sd_hsfr sd_hsne sd_hsmw _
sd_hssou sd_hswst sd_hsbr sd_hmob sd_pmnv sd_pmno sd_pmdel sd_mocmq sd_msondq sd_fsncom sd_fspcom _
sd_fspin sd_fspcap sd_fsput sd_fsdxp sd_fspxe sd_exrsw sd_exrjan sd_exruk sd_exrcan sd_fygm3 sd_fygm6 sd_fygt1 _
sd_fygt5 sd_fygt10 sd_fyaaac sd_fybaac sd_sfygm3 sd_sfygm6 sd_sfygt1 sd_sfygt5 sd_sfygt10 sd_sfyaaac sd_sfybaac _
sd_fm1 sd_fm2 sd_fm3 sd_fm2dq sd_fmfba sd_fmrra sd_fmrnba sd_fclnq sd_fclbmc sd_ccinrv sd_pmcp sd_pwfsa sd_pwfcsa _ 
sd_pwimsa sd_pwcmsa sd_psm99q sd_punew sd_pu83 sd_pu84 sd_pu85 sd_puc sd_pucd sd_pus sd_puxf sd_puxhs _
sd_puxm sd_lehcc sd_lehm sd_hhsntn

x_series_f3.makepcomp(eigval=x_eigval_f3,eigvec=x_eigvec_f3) chat_f1 chat_f2 chat_f3 chat_f4 chat_f5 'These components include all the series. The eigvectors are stored in the matrix all_eigvec

'The second step is to extract 5 principal components from the 70 X "slow moving" X variables.  The variables are ordered from slowest to fastest "moving"...

Group nxslow sd_ip sd_lhur sd_punew sd_ipp sd_ipf sd_ipc sd_ipcd sd_ipcn sd_ipe sd_ipi sd_ipm sd_ipmd sd_ipmnd sd_ipmfg _
sd_ipd sd_ipn sd_ipmin sd_iput sd_ipxmca sd_pmi sd_pmp sd_gmpyq sd_gmyxpq sd_lhel sd_lhelx sd_lhem sd_lhnag _
sd_lhu680 sd_lhu5 sd_lhu14 sd_lhu15 sd_lhu26 sd_lpnag sd_lp sd_lpgd sd_lpmi sd_lpcc sd_lpem sd_lped _
sd_lpen sd_lpsp sd_lptu sd_lpt sd_lpfr sd_lps sd_lpgov sd_lphrm sd_lpmosa sd_pmemp sd_gmcq sd_gmcdq sd_gmcnq _
sd_gmcsq sd_gmcanq sd_pwfsa sd_pwfcsa sd_pwimsa sd_pwcmsa sd_psm99q sd_pu83 sd_pu84 sd_pu85 sd_puc sd_pucd _
sd_pus sd_puxf sd_puxhs sd_puxm sd_lehcc sd_lehm 

'Extract the slow principal components

nxslow.makepcomp(eigval=slow_eigval_f3,eigvec=slow_eigvec_f3) slow_f1 slow_f2 slow_f3  slow_f4 slow_f5 'Data * slow_eigvec

'These principal components need to be regressed on fyff (R in the paper, see page 405, top paragraph) 

equation f1star.ls chat_f1 sd_fyff slow_f1 slow_f2 slow_f3 
genr f1 = chat_f1 - (f1star.@coefs(1)*sd_fyff)

equation f2star.ls chat_f2 sd_fyff slow_f1 slow_f2 slow_f3 
genr f2 = chat_f2 - (f2star.@coefs(1)*sd_fyff)

equation f3star.ls chat_f3  sd_fyff slow_f1 slow_f2 slow_f3
genr f3 = chat_f3 - (f3star.@coefs(1)*sd_fyff)

'Run the first VAR...baseline with just fyff

delete(noerr) imp_* var_0* var_1* 

var var_0.ls(noconst) 1 !var_lags sd_ip sd_punew sd_fyff
var_0.impulse(!imp_periods,imp=chol,matbys=imp_0) sd_ip sd_punew sd_fyff  @ sd_fyff

matrix imp_0 = (imp_0/imp_0(1,9))*!dampener

vector ip_1 = @columnextract(imp_0, 9)
smpl @last-!imp_periods+1 @last 
mtos(ip_1, f0_impulse_fyff)  'Convert the first column of imp_1 to a series...
line f0_impulse_fyff

vector ip_1 = @columnextract(imp_0, 8)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f0_impulse_punew) 'Convert the second column of imp_1 to a series...
line f0_impulse_punew
genr f0_impulse_punew = exp(@cumsum(f0_impulse_punew))-1

vector ip_1 = @columnextract(imp_0, 7)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f0_impulse_ip) 'Convert the third column of imp_1 to a series...
line f0_impulse_ip
genr f0_impulse_ip = exp(@cumsum(f0_impulse_ip))-1

smpl @all

var f3_var_1.ls(noconst) 1 !var_lags  f1 f2 f3 sd_fyff @ c

'Compute the impulse response matrix with respect to all the variables in the VAR

delete(noerr) f3_imp_1

f3_var_1.impulse(!imp_periods,imp=chol,matbys=f3_imp_1) f1 f2 f3 sd_fyff @ sd_fyff

matrix f3_imp_1 = (f3_imp_1/f3_imp_1(1,16))*!dampener  '(This number scales the impulse response functions so that the shock is 25 basis points)

show f3_imp_1

delete f3_impulse_*

vector ip_1 = @columnextract(f3_imp_1, 13)
smpl @last-!imp_periods+1 @last 
mtos(ip_1, f3_impulse_f1)  'Convert the first column of imp_1 to a series...
line f3_impulse_f1

vector ip_1 = @columnextract(f3_imp_1, 14)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f3_impulse_f2) 'Convert the second column of imp_1 to a series...
line f3_impulse_f2

vector ip_1 = @columnextract(f3_imp_1, 15)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f3_impulse_f3) 'Convert the third column of imp_1 to a series...
line f3_impulse_f3

vector ip_1 = @columnextract(f3_imp_1, 16)
smpl @last-!imp_periods+1 @last 
mtos(ip_1, f3_impulse_fyff) 'Convert the fourth column of imp_1 to a series...
line f3_impulse_fyff

smpl @first @last-!imp_periods

genr f3_impulse_f1 = 0
genr f3_impulse_f2 = 0
genr f3_impulse_f3 = 0
genr f3_impulse_fyff = 0

smpl @all

'Now regress the X variables of interest on the original factors and sd_fyff to derive the impulse response function of the X variables of interest

for %s ip punew fm2 fmfba exrjan gmcq gmcdq gmcnq lehcc
	smpl @all
	equation x_impulse.ls sd_{%s} f1 f2 f3 sd_fyff
	genr f3_impulse_{%s} = x_impulse.@coefs(1)*f3_impulse_f1 + x_impulse.@coefs(2)*f3_impulse_f2 + x_impulse.@coefs(3)*f3_impulse_f3+x_impulse.@coefs(4)*f3_impulse_fyff
     smpl @last-!imp_periods+1 @last
    	genr f3_impulse_{%s} = exp(@cumsum(f3_impulse_{%s}))-1
 next 

for %s hsfr
     smpl @all
	equation x_impulse.ls sd_{%s} f1 f2 f3 sd_fyff
	genr f3_impulse_{%s} = x_impulse.@coefs(1)*f3_impulse_f1 + x_impulse.@coefs(2)*f3_impulse_f2 + x_impulse.@coefs(3)*f3_impulse_f3+x_impulse.@coefs(4)*f3_impulse_fyff
	smpl @last-!imp_periods+1 @last
    	genr f3_impulse_{%s} = exp(f3_impulse_{%s})-1
next 

for %s fygm3 fygt5 pmcp ipxmca lhur pmemp pmno fsdxp hhsntn 
	smpl @all
	equation x_impulse.ls sd_{%s} f1 f2 f3 sd_fyff
	genr f3_impulse_{%s} = x_impulse.@coefs(1)*f3_impulse_f1 + x_impulse.@coefs(2)*f3_impulse_f2 + x_impulse.@coefs(3)*f3_impulse_f3+x_impulse.@coefs(4)*f3_impulse_fyff
 next  

smpl @last-!imp_periods+1 @last

delete(noerr) bbe_figure_2

graph bbe_figure_2.line(m) f3_impulse_fyff f3_impulse_ip f3_impulse_punew f3_impulse_fygm3 f3_impulse_fygt5 f3_impulse_fmfba f3_impulse_fm2 f3_impulse_exrjan f3_impulse_pmcp f3_impulse_ipxmca f3_impulse_gmcq f3_impulse_gmcdq f3_impulse_gmcnq f3_impulse_lhur f3_impulse_pmemp f3_impulse_lehcc f3_impulse_hsfr f3_impulse_pmno f3_impulse_fsdxp f3_impulse_hhsntn
bbe_figure_2.axis(l) zeroline -minor font(12)
bbe_figure_2.addtext(t,font(,30,b)) Figure II: Impulse Response Functions from FAVAR using Three Factors (Principal Components)

delete(noerr) graph_f1

graph graph_f1.line(m) f1_impulse_fyff f1_impulse_ip f1_impulse_punew f1_impulse_fygm3 f1_impulse_fygt5 f1_impulse_fmfba f1_impulse_fm2 f1_impulse_exrjan f1_impulse_pmcp f1_impulse_ipxmca f1_impulse_gmcq f1_impulse_gmcdq f1_impulse_gmcnq f1_impulse_lhur f1_impulse_pmemp f1_impulse_lehcc f1_impulse_hsfr f1_impulse_pmno f1_impulse_fsdxp f1_impulse_hhsntn

graph_f1.axis(l) zeroline -minor font(12)
graph_f1.addtext(t,font(,30,b)) Figure I: Impulse Response Functions from FAVAR using one Factor, CPI and IP (Principal Components)

delete(noerr) graph_f0_*

graph graph_f0_cpi.line f0_impulse_punew f1_impulse_punew f3_impulse_punew
graph_f0_cpi.axis(l) zeroline -minor font(12)
graph graph_f0_ip.line f0_impulse_ip f1_impulse_ip f3_impulse_ip
graph_f0_ip.axis(l) zeroline -minor font(12)
graph graph_f0_fyff.line f0_impulse_fyff f1_impulse_fyff f3_impulse_fyff
graph_f0_fyff.axis(l) zeroline -minor font(12)

delete(noerr) bbe_figure_1*

graph bbe_figure_1.merge graph_f0_fyff graph_f0_ip graph_f0_cpi
bbe_figure_1.addtext(t,font(,30,b)) Figure I: Impulse Response Functions from FAVAR using one Factor, CPI and IP (Principal Components)
