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

group x_series_f1 sd_ipp sd_ipf sd_ipc sd_ipcd sd_ipcn sd_ipe sd_ipi sd_ipm sd_ipmd sd_ipmnd sd_ipmfg sd_ipd sd_ipn sd_ipmin sd_iput _
sd_ipxmca sd_pmi sd_pmp sd_gmpyq sd_gmyxpq sd_lhel sd_lhelx sd_lhem sd_lhnag sd_lhur sd_lhu680 sd_lhu5 sd_lhu14 _
sd_lhu15 sd_lhu26 sd_lpnag sd_lp sd_lpgd sd_lpmi sd_lpcc sd_lpem sd_lped sd_lpen sd_lpsp sd_lptu sd_lpt sd_lpfr sd_lps sd_lpgov _
sd_lphrm sd_lpmosa sd_pmemp sd_gmcq sd_gmcdq sd_gmcnq sd_gmcsq sd_gmcanq sd_hsfr sd_hsne sd_hsmw _
sd_hssou sd_hswst sd_hsbr sd_hmob sd_pmnv sd_pmno sd_pmdel sd_mocmq sd_msondq sd_fsncom sd_fspcom _
sd_fspin sd_fspcap sd_fsput sd_fsdxp sd_fspxe sd_exrsw sd_exrjan sd_exruk sd_exrcan sd_fygm3 sd_fygm6 sd_fygt1 _
sd_fygt5 sd_fygt10 sd_fyaaac sd_fybaac sd_sfygm3 sd_sfygm6 sd_sfygt1 sd_sfygt5 sd_sfygt10 sd_sfyaaac sd_sfybaac _
sd_fm1 sd_fm2 sd_fm3 sd_fm2dq sd_fmfba sd_fmrra sd_fmrnba sd_fclnq sd_fclbmc sd_ccinrv sd_pmcp sd_pwfsa sd_pwfcsa _ 
sd_pwimsa sd_pwcmsa sd_psm99q sd_pu83 sd_pu84 sd_pu85 sd_puc sd_pucd sd_pus sd_puxf sd_puxhs _
sd_puxm sd_lehcc sd_lehm sd_hhsntn

x_series_f1.makepcomp(eigval=x_eigval_f1,eigvec=x_eigvec_f1) chat_f1  'These components include all the series except fyff, ip and punew. The eigvectors are stored in the matrix x_eigvec

'The second step is to extract 1 principal components from the 68 X "slow moving" X variables.  The variables are ordered from slowest to fastest "moving"...

Group nxslow_f1 sd_lhur sd_ipp sd_ipf sd_ipc sd_ipcd sd_ipcn sd_ipe sd_ipi sd_ipm sd_ipmd sd_ipmnd sd_ipmfg _
sd_ipd sd_ipn sd_ipmin sd_iput sd_ipxmca sd_pmi sd_pmp sd_gmpyq sd_gmyxpq sd_lhel sd_lhelx sd_lhem sd_lhnag _
sd_lhu680 sd_lhu5 sd_lhu14 sd_lhu15 sd_lhu26 sd_lpnag sd_lp sd_lpgd sd_lpmi sd_lpcc sd_lpem sd_lped _
sd_lpen sd_lpsp sd_lptu sd_lpt sd_lpfr sd_lps sd_lpgov sd_lphrm sd_lpmosa sd_pmemp sd_gmcq sd_gmcdq sd_gmcnq _
sd_gmcsq sd_gmcanq sd_pwfsa sd_pwfcsa sd_pwimsa sd_pwcmsa sd_psm99q sd_pu83 sd_pu84 sd_pu85 sd_puc sd_pucd _
sd_pus sd_puxf sd_puxhs sd_puxm sd_lehcc sd_lehm 

'Extract the slow principal components

nxslow_f1.makepcomp(eigval=slow_eigval_f1,eigvec=slow_eigvec_f1) slow_f1

'These principal components need to be regressed on fyff (R in the paper, see page 405, top paragraph) 

equation f1star.ls chat_f1 sd_fyff slow_f1 
genr f1 = chat_f1 - (f1star.@coefs(1)*sd_fyff)

'Run the first VAR...baseline with just fyff

delete f1_impulse_* f1_var_0* f1_var_1* f1_imp_*

smpl @all

var f1_var_1.ls(noconst) 1 !var_lags f1 sd_ip sd_punew sd_fyff

'Compute the impulse response matrix with respect to all the variables in the VAR

f1_var_1.impulse(!imp_periods,imp=chol,matbys=f1_imp_1) f1 sd_ip sd_punew sd_fyff @ sd_fyff

matrix f1_imp_1 = (f1_imp_1/f1_imp_1(1,16))*!dampener  '(This number scales the impulse response functions so that the shock is 25 basis points)

delete(noerr) f1_impulse_*

vector ip_1 = @columnextract(f1_imp_1, 13)
smpl @last-!imp_periods+1 @last 
mtos(ip_1, f1_impulse_f1)  'Convert the first column of imp_1 to a series...
line f1_impulse_f1

vector ip_1 = @columnextract(f1_imp_1, 14)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f1_impulse_ip) 'Convert the second column of imp_1 to a series...
line f1_impulse_ip

vector ip_1 = @columnextract(f1_imp_1, 15)
smpl @last-!imp_periods+1 @last
mtos(ip_1, f1_impulse_punew) 'Convert the third column of imp_1 to a series...
line f1_impulse_punew

vector ip_1 = @columnextract(f1_imp_1, 16)
smpl @last-!imp_periods+1 @last 
mtos(ip_1, f1_impulse_fyff) 'Convert the fourth column of imp_1 to a series...
line f1_impulse_fyff

smpl @first @last-!imp_periods

genr f1_impulse_fyff = 0
genr f1_impulse_punew = 0
genr f1_impulse_ip= 0
genr f1_impulse_fyff = 0

smpl @all

'Now regress the X variables of interest on the original factors and sd_fyff to derive the impulse response function of the X variables of interest

for %s fm2 fmfba exrjan gmcq gmcdq gmcnq lehcc
	smpl @all
	equation x_impulse.ls sd_{%s} f1 sd_fyff
	genr f1_impulse_{%s} = x_impulse.@coefs(1)*f1_impulse_f1+x_impulse.@coefs(2)*f1_impulse_fyff
     smpl @last-!imp_periods+1 @last
    	genr f1_impulse_{%s} = exp(@cumsum(f1_impulse_{%s}))-1
 next 

for %s ip punew
     smpl @last-!imp_periods+1 @last
    	genr f1_impulse_{%s} = exp(@cumsum(f1_impulse_{%s}))-1
next

for %s hsfr
     smpl @all
	equation x_impulse.ls sd_{%s} f1 sd_fyff
	genr f1_impulse_{%s} = x_impulse.@coefs(1)*f1_impulse_f1 + x_impulse.@coefs(2)*f1_impulse_fyff
	smpl @last-!imp_periods+1 @last
    	genr f1_impulse_{%s} = exp(f1_impulse_{%s})-1
next 

for %s fygm3 fygt5 pmcp ipxmca lhur pmemp pmno fsdxp hhsntn 
	smpl @all
	equation x_impulse.ls sd_{%s} f1 sd_fyff
	genr f1_impulse_{%s} = x_impulse.@coefs(1)*f1_impulse_f1+x_impulse.@coefs(2)*f1_impulse_fyff
next


