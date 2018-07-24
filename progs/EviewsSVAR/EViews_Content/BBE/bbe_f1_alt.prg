pageselect starthere 'Always work with the pagefile called "starthere"

%cs = @pagesmpl  'Make a copy of the current sample to restore it later...

smpl 1959:01 2001:7 'Set the working sample...

scalar nvars = 120 'Total number of variables in the dataset, including the Federal Funds Rate

scalar nslow = 70  'Number of variables in the slow moving dataset...

!imp_periods = 48 'Length of the impulse response functions...

%match_bbe = "y"  'Set to "y" to match the scale of the impulse response functions in Bernanke et. al.

delete(noerr) x_slow* x_fast* x_all* 'Start from a fresh slate by removing existing versions of the PC's

%ss = "fyff"  'Name the variables to be excluded from the PC list here, separating each excluded variable with an underscore ("fyff_ipe")...

group x_all  'Given the current exclusion list, this group (x_all) will eventually contain 119 variables at the end of the run; basically all the data except sd_fyff

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

     if @instr(%ss,%s) = 0 then 
        	x_all.add sd_{%s}
	endif

next

!dampener = 0.25/@stdev(fyff)  '~8 basis points

'The first step is to extract K principal components from the 70 X "slow moving" X variables.  The variables are ordered from slowest to fastest moving...

Group x_slow sd_ip sd_lhur sd_ipp sd_ipf sd_ipc sd_ipcd sd_ipcn sd_ipe sd_ipi sd_ipm sd_ipmd sd_ipmnd sd_ipmfg _
sd_ipd sd_ipn sd_ipmin sd_iput sd_ipxmca sd_pmi sd_pmp sd_gmpyq sd_gmyxpq sd_lhel sd_lhelx sd_lhem sd_lhnag _
sd_lhu680 sd_lhu5 sd_lhu14 sd_lhu15 sd_lhu26 sd_lpnag sd_lp sd_lpgd sd_lpmi sd_lpcc sd_lpem sd_lped _
sd_lpen sd_lpsp sd_lptu sd_lpt sd_lpfr sd_lps sd_lpgov sd_lphrm sd_lpmosa sd_pmemp sd_gmcq sd_gmcdq sd_gmcnq _
sd_gmcsq sd_gmcanq sd_pwfsa sd_pwfcsa sd_pwimsa sd_pwcmsa sd_psm99q sd_pu83 sd_pu84 sd_pu85 sd_puc sd_pucd _
sd_pus sd_puxf sd_puxhs sd_puxm sd_lehcc sd_lehm  sd_punew

'These are the fast moving variables...

group x_fast sd_hsfr	sd_hsne	sd_hsmw	sd_hssou	sd_hswst	sd_hsbr	sd_hmob	sd_pmnv	sd_pmno	sd_pmdel	sd_mocmq	sd_msondq	sd_fsncom	sd_fspcom	sd_fspin	sd_fspcap	sd_fsput	sd_fsdxp	sd_fspxe	sd_exrsw	sd_exrjan	sd_exruk	sd_exrcan	sd_fygm3	sd_fygm6	sd_fygt1	sd_fygt5	sd_fygt10	sd_fyaaac	sd_fybaac	sd_sfygm3	sd_sfygm6	sd_sfygt1	sd_sfygt5	sd_sfygt10	sd_sfyaaac	sd_sfybaac	sd_fm1	sd_fm2	sd_fm3	sd_fm2dq	sd_fmfba	sd_fmrra	sd_fmrnba	sd_fclnq	sd_fclbmc	sd_ccinrv	sd_pmcp	sd_hhsntn

'Extract the slow principal components. Use !k_ps below to change the number of PCs used...

!k_ps = 3

%slow_list = "slow_f1 "

for !j = 2 to !k_ps
	%slow_list = %slow_list + "slow_f"  + @str({!j}) + " "
next

%slow_list = @trim(%slow_list)

x_slow.makepcomp(eigval=slow_eigval_f1,eigvec=slow_eigvec_f1) {%slow_list}

'Now estimate an SVAR involving the slow moving PC's and the federal funds rate...The PC(s) component is a standard VAR...

!n_lags = 13 'Number of lags used to estimate the SVAR

var svar_pc.ls 1 !n_lags {%slow_list} sd_fyff

delete(noerr) pat_a* pat_b* f1_imp  'Remove existing pattern matrices.  Start from a clean slate...

'Create the "A" pattern matrix

matrix pat_a = @identity(!k_ps+1)
for !jj = 1 to !k_ps
	pat_a(!k_ps+1,!jj) = na  'Basically the identity matrix except for the last row, which has unknown parameters
next

'Create the "B" pattern matrix

matrix pat_b = @identity(!k_ps+1)
for !jj = 1 to !k_ps+1
	pat_b(!jj,!jj) = na 'Diagonal matrix ... 0's on the off-diagonals; variances on the diagonal terms.
next

'Estimate the SVAR given the current pat_a and pat_b matrix

svar_pc.svar(rtype=patsr,namea=pat_a,nameb=pat_b,f0=u,nostop,conv=1e-10,maxiter=2500)

if svar_pc.@svarcvgtype <> 0 then 'This line tests whether the SVAR converged.  It ends the run if not...typically due to bad starting values.

	show svar_pc

	@uiprompt("Structural VAR: structural VAR failed to converge!")

	stop

endif

!n = !imp_periods

'Kill the previous impulse response function set...namely, f1_imp

delete(noerr) f1_imp

svar_pc.impulse(!imp_periods,matbys=f1_imp,se=a,imp=struct) {%slow_list} sd_fyff

'Given the use of the matbys option to order the impulse responses, the first column of "f1_imp" is the response of the first variable to the first shock, the second column is the response of the second variable to the first shock, and so on. The response and shock orderings correspond to the ordering of variables in the VAR.

'We are therefore interested in the last (!k_ps+1) columns, given there is only one endogenous variable...fyff

'The impulse response functions are saved in a matrix called "rel_impulses"

'show f1_imp

!jk = (!k_ps+1)^2 'Total number of impulses = number of variables^2

matrix rel_impulses = @subextract(f1_imp,1,!jk-(!k_ps),!n,!jk)

'show rel_impulses

rel_impulses = (rel_impulses/f1_imp(1,!jk))*!dampener  '(This scales the impulse response functions so that the R shock is 25 basis points)

'Copy the impulse response function set to matrices (vectors, actually) r1, r2, r3, ...,

for !j = 1 to !k_ps+1
	matrix r{!j} = @subextract(rel_impulses,1,!j,!n,!j)
next

delete(noerr) s_response_* r_response_* a_response f_response_*

'First, estimate the response of the slow moving variables to the principal components...

%slow_names = x_slow.@members

for %s {%slow_names}
	smpl @all
	equation x_impulse_{%s}.ls {%s} {%slow_list} 'Regress the slow moving variable just on the PCs
	matrix s_response_{%s} = r1*x_impulse_{%s}.@coefs(1)
	for !j = 2 to !k_ps
	 	matrix s_response_{%s} = s_response_{%s} + (r{!j}*x_impulse_{%s}.@coefs(!j))
	next
     matrix a_response_{%s} = s_response_{%s} 'Store the final result in variables that begin with a_..."
next 

'Second, estimate the response of the fast moving variables to the principal components...

%fast_names = x_fast.@members

for %s {%fast_names}
	smpl @all
	equation x_impulse_{%s}.ls {%s} {%slow_list} sd_fyff 'Regress the fast moving variable on the principal components AND sd_fyff
	matrix f_response_{%s} = r1*x_impulse_{%s}.@coefs(1)
	for !j = 2 to !k_ps+1
		matrix f_response_{%s} = f_response_{%s} + (r{!j}*x_impulse_{%s}.@coefs(!j))
	next
	matrix a_response_{%s} = f_response_{%s}
next

delete(noerr) p1* gr_* gr_all*  'Clean up the old graphs/variables...

'Now build the comparative graphs...

 %grlist = ""

'Note: the comparison will compare like with like; the impulse response functions in the workfile (calculated by the original replication program) in Bernanke are in terms of "exp(@cumsum(p1))-1". Thus the use of dlog(f1_impulse_{%s}+1) below 

 for %s fm2 fmfba exrjan gmcq gmcdq gmcnq lehcc ip punew
	series p1
	smpl @last-!imp_periods+1 @last
	mtos(a_response_sd_{%s}, p1)
     if %match_bbe = "y" then
   		genr pagan_f_impulse_{%s} = exp(@cumsum(p1))-1 '-- Original expression taken from Bernanke, et. al. 
     		graph gr_{%s} pagan_f_impulse_{%s}  f{!k_ps}_impulse_{%s}
     		%tt = "Response of " + @upper(%s) + " to an interest rate shock"
	else
    		genr pagan_f_impulse_{%s} = 100*p1  
     		graph gr_{%s} pagan_f_impulse_{%s}  100*dlog(f{!k_ps}_impulse_{%s}+1) 'Need to unwind the Bernanke expression...
     		%tt = "Response (x 100) of " + @upper(%s) + " to an interest rate shock"
	endif
     '%tt = "Response (x 100) of " + @upper(%s) + " to a 25 basis point increase in R"
     gr_{%s}.addtext(t, font("Times",16)) {%tt}
	gr_{%s}.setelem(2) lpat(dash6)
	gr_{%s}.options -color
     %grlist = %grlist + " gr_"+%s
     delete(noerr) p1
  next

 graph gr_first.merge {%grlist}
 gr_first.addtext(t, font("Times",40)) Pagan versus Bernake et al
 show gr_first

 @uiprompt("Press OK to see the next Pagan versus Benanke graph!","O")

 for %s hsfr
	series p1
	smpl @last-!imp_periods+1 @last
	mtos(a_response_sd_{%s}, p1)
     if %match_bbe = "y" then
   		genr pagan_f_impulse_{%s} = exp(p1)-1 '-- Original expression taken from Bernanke, et. al.; necessary because they want to compare impulse to the level of the original series.  Note that the raw series was expressed in log() form 
     		graph gr_{%s} pagan_f_impulse_{%s}  f{!k_ps}_impulse_{%s}
     		%tt = "Response of " + @upper(%s) + " to an interest rate shock"
	else
    		genr pagan_f_impulse_{%s} = 100*p1
     		graph gr_{%s} pagan_f_impulse_{%s}  100*log(f{!k_ps}_impulse_{%s}+1) 'Unwind BBE transform...
     		%tt = "Response (x 100) of " + @upper(%s) + " to an interest rate shock"
	endif
	gr_{%s}.setelem(2) lpat(dash6)
	gr_{%s}.options -color
     gr_{%s}.addtext(t, font("Times",16)) {%tt}
     delete(noerr) p1
 next 

 show gr_hsfr
 @uiprompt("Press OK to see the next Pagan versus Bernanke graph!","O")

 %grlist = ""

 for %s fygm3 fygt5 pmcp ipxmca lhur pmemp pmno fsdxp hhsntn 
	series p1
	smpl @last-!imp_periods+1 @last
	mtos(a_response_sd_{%s}, p1)
     if %match_bbe = "y" then
   		genr pagan_f_impulse_{%s} = p1 
     		graph gr_{%s} pagan_f_impulse_{%s}  f{!k_ps}_impulse_{%s} 'Impulse response functions are directly comparable for these variables. 
          %tt = "Response of " + @upper(%s) + " to an interest rate shock"
	else
    		genr pagan_f_impulse_{%s} = 100*p1
     		graph gr_{%s} pagan_f_impulse_{%s}  100*f{!k_ps}_impulse_{%s} 'Impulse response functions are directly comparable...
          %tt = "Response (x 100) of " + @upper(%s) + " to an interest rate shock"
	endif
	gr_{%s}.setelem(2) lpat(dash6)
	gr_{%s}.options -color
     gr_{%s}.addtext(t, font("Times",16)) {%tt}
     %grlist = %grlist + " gr_"+%s
     delete(noerr) p1
 next

 graph gr_second.merge {%grlist}
 gr_second.addtext(t, font("Times",40)) Pagan versus Bernake et al
 show gr_second

'Lastly, restrict the graphs to those selected for the OPR paper...

for %s lehcc GMCDQ PUNEW EXRJAN IP HSFR

	gr_{%s}.name(1) OPR

	gr_{%s}.name(2) BBE

	gr_{%s}.axis(b) angle(0)

	gr_{%s}.axis(b) font("Times",6)

	for !kk = 1 to !n

		gr_{%s}.setobslabel(!kk) " "
 
		if @mod(!kk , 5) = 0 then 

		   gr_{%s}.setobslabel(!kk) {!kk}

		endif

	next

	gr_{%s}.setelem(2) lpat(dash6)

	gr_{%s}.options -color

	gr_{%s}.setobslabel(1) 1

	gr_{%s}.setobslabel(!n) {!n}

next

'gr_LEHCC.addtext(t,font("Times",12)) AVG HR EARNINGS OF CONSTR WKRS: CONSTRUCTION ($,SA)

'gr_GMCDQ.addtext(t,font("Times",12)) PERSONAL CONS. EXP (CHAINED) – TOT. DUR. (BIL 96$,SAAR)

'gr_PUNEW.addtext(t,font("Times",12)) CPI-U: ALL ITEMS (82-84=100,SA)

'gr_EXRJAN.addtext(t,font("Times",12)) AVG HR EARNINGS OF CONSTR WKRS: CONSTRUCTION ($,SA)

'gr_IP.addtext(t,font("Times",12)) INDUSTRIAL PRODUCTION: TOTAL INDEX (1992=100,SA)

'gr_HSFR.addtext(t,font("Times",12)) HOUSING STARTS:NONFARM(1947-58);TOT.(1959-) (THOUS.,SA)

graph gr_last.merge gr_lehcc gr_GMCDQ gr_PUNEW gr_EXRJAN gr_IP gr_HSFR

show gr_last

