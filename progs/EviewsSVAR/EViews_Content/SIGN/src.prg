smpl @all

pageselect chomoreno

delete(noerr) impulse_set_* counter* gp_* gr_* i*_* svar_impulse* sign_mat* inst* e1 e2 e3 a_* v* 

!nvars = 3         	'Number of variables in the base VAR

!n = 32 			'Length of a single impulse response function

!window = !n		'Window: criterion: first eight values of the impulse response vector must have the correct sign; not used for this implementation

!j   = (!nvars^2) 	'Number of impulse response vectors in a single run

!niters =  5000000   'Total number of impulse response function sets to generate; number taken will typically be much smaller.

!number_required = 1000 'Maximum number required...

rndseed 123456789

for !s = 1 to !j
	delete(noerr) impulse_set_{!s}
	matrix(!n,!niters) impulse_set_{!s}=na  'This will ultimately store all the valid impuse response vectors
	%vnames_{!s} = "" 'Keeps track of the names of these vectors...
next

delete(noerr) sign_matrix ch

matrix(!nvars,!nvars) pattern_a = @identity(!nvars)
matrix(!nvars,!nvars) pattern_b = 0

scalar success = 0

var ch.ls 1 2 gap infl ff 'Create the VAR

smpl 1981q3 2000q1

ch.fit _f

delete(noerr) vf

group vf gap_f infl_f ff_f

stom(vf, m_vf)

'show m_vf

%smpl = ch.@smpl

scalar counter = 0  'This will keep track of the valid impulse count for each shock.

for !x = 1 to !niters 'Main iteration loop

	'Generate a new set of identifying parameters...

	pattern_a = @identity(!nvars)

	smpl {%smpl}

	!p1 = @runif(-1,1) 'Random parameters

	!p2 = @runif(-1,1)

	!p3 = @runif(-1,1)

	'Update the pattern matrix [created previously -- see pattern_a and pattern_b]

	pattern_a(1,2) = !p1/(1-@abs(!p1))  'Randomize the parameter

	pattern_a(1,3) = !p2/(1-@abs(!p2))  'Ditto

	pattern_a(2,3) = !p3/(1-@abs(!p3))  'Ditto

	'Estimate the remaining unknown parameters using IV so as to get unbiased estimates

	'Gap equation...

	genr y_gap = gap - pattern_a(1,2)*infl - pattern_a(1,3)*ff

	equation instr_gap.tsls y_gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c @ gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c

	genr e_gap = resid

	pattern_b(1,1) = instr_gap.@se

	'Inflation equation...

	equation instr_infl.tsls (infl + pattern_a(2,3)*ff) gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)  c @ e_gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c

	pattern_a(2,1) = -instr_infl.@coefs(1)

	genr e_infl = resid

	pattern_b(2,2) = instr_infl.@se^2

	'Interest rate equation

	equation instr_ff.tsls ff infl gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c @ e_gap e_infl gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c

	genr e_ff = resid

	pattern_a(3,2) = -instr_ff.@coefs(1)

	pattern_a(3,1) = -instr_ff.@coefs(2)

	pattern_b(3,3) = instr_ff.@se

	scalar x_3 = 1.0 - ((pattern_b(3,3)^2)/@vars(ff))

	matrix a_inv = @inverse(pattern_a)

	delete(noerr) ch_svar_imp* 'Remove the previous impulse response function set

	'freeze(mode=overwrite,impgraph) ch.impulse(!n,matbys=ch_svar_imp,se=a,imp=user,fname=shock_1) 

	'Sign restriction test ... c0 is the first row of ch_svar_imp

	'Check for a positive or negative demand shock on the basis of the inverse of the current pattern_a matrix

	!test_p = -1 'Assume no positive shock

	!test_n = -1 'Assume no negative shock

	!bad_draw = -1 '-1 is false

	!col_position = -1

	matrix(3,2) col_type = 0

	!pd_count = 0
	!nd_count = 0
	!pp_count = 0
	!np_count = 0
	!pr_count  = 0
	!nr_count = 0

	!t1 = 0
	!t2 = 0
	!t3 = 0

	!d_column = 0
	!p_column = 0
	!r_column = 0

	'Demand Shock

	'Notice the test bases the assessment just on a_inv and the signs of the coefficients in the matrix.  The reason is that the sign of the impulse response function at time t = 0 is implied.

	for !j = 1 to 3

		!test_p = ((a_inv(1,!j) >= 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) >= 0)) ' + Demand shock...

		if !test_p > 0 then 

			!pd_count = !pd_count + 1

			!d_column = !j

		endif

		!test_n = ((a_inv(1,!j) < 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) < 0))  ' - Demand shock...

		if !test_n > 0 then 

			!nd_count = !nd_count + 1

			!d_column = -!j  'Store the negative of the column index (to flag a negative demand shock); adjusted later.

		endif

	next

	if ((!pd_count+!nd_count) = 1) then 

		!t1 = 1 'This implies that there was a demand shock (positive or negative)

	endif

	'Inflation

	for !j = 1 to 3

		!test_p = ((a_inv(1,!j) < 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) >= 0))  '+Inflation shock

		if !test_p > 0 then 

			!pp_count = !pp_count + 1

			!p_column = !j

		endif

		!test_n = ((a_inv(1,!j) >= 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) < 0)) '- Inflation shock

		if !test_n > 0 then 

			!np_count = !np_count + 1

			!p_column = -!j

		endif

	next

	if ((!pp_count+!np_count) = 1) then 

		!t2 = 1 'This implies that there was an inflation shock (positive or negative)

	endif

	'Interest rate shock...

	for !j = 1 to 3

		!test_p = (a_inv(1,!j) < 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) >= 0) '+ interest rate shock

		if !test_p > 0 then 

			!pr_count = !pr_count + 1
	
			!r_column = !j

		endif

		!test_n = (a_inv(1,!j) >= 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) < 0) '- interest rate shock

		if !test_n > 0 then 

			!nr_count = !nr_count + 1

			!r_column = -!j

		endif

	next

	if ((!pr_count+!nr_count) = 1) then 

		!t3 = 1 'This implies an interest rate shock (+, -);

	endif

	if ((!t1+!t2 +!t3)>2.5) then 'At this point, we have one of each shock.  Re-create the A matrix so that they are ordered corrrectly.
	
		success = success + 1

'		Re-order the columns

		matrix final_a = a_inv

		matrix b_shocks = pattern_b

		vector v1 = @subextract(a_inv,1,@abs(!d_column),!nvars,@abs(!d_column))

		if !d_column < 0 then
		 	v1 = -v1
		endif

		colplace(final_a,v1,1)

		b_shocks(1,1) = pattern_b(@abs(!d_column),@abs(!d_column))

		vector v2 = @subextract(a_inv,1,@abs(!p_column),!nvars,@abs(!p_column))

		if !p_column < 0 then
			v2 = -v2
		endif

		colplace(final_a,v2,2)

		b_shocks(2,2) = pattern_b(@abs(!p_column),@abs(!p_column))

		vector v3 = @subextract(a_inv,1,@abs(!r_column),!nvars,@abs(!r_column))

		if !r_column < 0 then
			v3 = -v3
		endif

		colplace(final_a,v3,3)

		b_shocks(3,3) = pattern_b(@abs(!r_column),@abs(!r_column))

		matrix shock_1 = final_a '*b_shocks

'		Calculate the new impulse response function...

'		Impulse reponse functions are ordered as gap w.r.t gap, infl w.r.t gap, ff with respect to gap
'																gap w.r.t infl,   infl w.r.t. infl; ff w.r.t. infl
'																gap w.r.t ff, infl w.r.t ff, ff w.r.t ff	
'							

		freeze(mode=overwrite,impgraph) ch.impulse(!n,smat=ch_svar_imp,se=a,imp=user,fname=shock_1) gap infl ff @ gap infl ff @ gap infl ff
		
		counter = counter+1  'Update the valid counter index 
	
		%tt = @str(counter)

		%tt ="Success count = " + @str(%tt) + " models; " + @str((counter/!x)*100) + " percent. Counter = " + @str(!x) + "." 

		statusline {%tt}

		for !v = 1 to !nvars^2

			delete(noerr) p1

			smpl @last - !n + 1 @last

			!t = counter

			vector p1 = @columnextract(ch_svar_imp,!v)

			colplace(impulse_set_{!v},p1,!t)

			series i{!v}_{!t} = 0 'In addition to the impulse_set store, create a series for the valid impulse response function

			mtos(p1,i{!v}_{!t}) 'Save the valid impulse to a series (for graphing purposes)
                                                                       		
			%vnames_{!v} = %vnames_{!v} + " " + "i"+@str({!v})+"_"+@str({!t}) 'Keep track of the valid names...

		next

		'At this point, let's see what the actual forecasts look like..

		smpl 1981q1 2000q1

		matrix ahat_c = @inverse(final_a)

		genr gap_c = gap - ahat_c(1,2)*infl - ahat_c(1,3)*ff

		equation ols_c1.ls gap_c c gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)

		ols_c1.fit gap_f

		genr gap_f = gap_f + ahat_c(1,2)*infl + ahat_c(1,3)*ff

		smpl 1981q1 2000q1

		genr infl_c = infl - ahat_c(2,1)*gap - ahat_c(2,3)*ff

		equation ols_c2.ls infl_c c gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)

		ols_c2.fit infl_f

		genr infl_c = infl_f + ahat_c(2,1)*gap + ahat_c(2,3)*ff

		smpl 1981q1 2000q1

		genr ff_c = ff - ahat_c(3,1)*gap - ahat_c(3,2)*infl

		equation ols_c3.ls ff_c c gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)

		ols_c3.fit ff_f

		genr ff_c = ff_f + ahat_c(3,1)*gap + ahat_c(3,2)*ff

		'show gap_c infl_c ff_c

	endif

	!x_save = !x

	if (counter = !number_required) then 
		exitloop 'Reached the required amount...
	endif

next 'End of iteration loop

'Calculate the MT response...

'Retrieve first line of impulse response function set

smpl @last - !n + 1 @last

for !v = 1 to !nvars^2

	impulse_set_{!v} = @transpose(@subextract(impulse_set_{!v},1,1,!n,counter))

next

delete(noerr) vv, b_star_store

matrix(counter, !nvars^2) b_star_store = 0

for !v = 1 to !nvars^2

	colplace(b_star_store,@columnextract(impulse_set_{!v},1),!v) 

next

delete(noerr) med

matrix(!nvars^2,1) med = 0

matrix(!nvars^2,1) stdev = 0

matrix(counter,!nvars^2) pk = 0

for !ss = 1 to !nvars^2 'Across all the separate impulse response functions...

	matrix vv = @columnextract(b_star_store,!ss)

	med(!ss, 1) = @median(vv)

	stdev(!ss,1) = @stdev(vv)

	for !cc = 1 to counter

		pk(!cc,!ss) = (impulse_set_{!ss}(!cc,!ss)-med(!ss,1))/stdev(!ss,1) 

	next

next

delete(noerr) xmin

matrix(counter,1) xmin = 0

!min_index = 0

!oo = 10^12

for !cc = 1 to counter

	matrix r_pk = @rowextract(pk,!cc)

	!oox = r_pk*@transpose(r_pk)

	xmin(!cc,1) = !oox

	if xmin(!cc,1) < !oo then 

		!min_index = !cc
		
		!oo = xmin(!cc,1)

	endif

next

@uiprompt("Minimum index = " + @str(!min_index))

'Check for uniqueness....

!nn = 0

for !v = 1 to counter

	if xmin(!v,1) = !oo then 

		!nn = !nn + 1

	endif

next

if (!nn > 1) then 

	@uiprompt("Warning: MT Criterion: uniqueness test failed.")

endif

if (!nn < 1) then 

	@uiprompt("Warning: MT Criterion: logic error.")

endif

'Now that we have the minimum index, save the impulse response functions in the appropriate location

for !v = 1 to !nvars^2
	impulse_set_{!v} = @transpose(impulse_set_{!v})
	matrix c_temp = @columnextract(impulse_set_{!v},!min_index)
	delete(noerr) median_imp_s_{!v}
	series median_imp_s_{!v} = 0
	mtos(c_temp, median_imp_s_{!v})
next

for !v = 1 to !nvars^2
	delete(noerr) med_gr_s_{!v}
	graph med_gr_s_{!v}.line median_imp_s_{!v}
	median_imp_s_{!v}.label SRC
	med_gr_s_{!v}.options -legend
	med_gr_s_{!v}.axis(l) zeroline
next

med_gr_s_1.addtext(t, font("Times",20)) Output Gap Response to a Demand Shock
med_gr_s_2.addtext(t, font("Times",20)) Inflation Response to a Demand Shock
med_gr_s_3.addtext(t, font("Times",20)) Interest Rate Response to a Demand Shock
med_gr_s_4.addtext(t, font("Times",20)) Output Gap Response to a Cost-Push Shock
med_gr_s_5.addtext(t, font("Times",20)) Inflation Response to a Cost-Push Shock
med_gr_s_6.addtext(t, font("Times",20)) Interest Rate Response to a Cost-Push Shock
med_gr_s_7.addtext(t, font("Times",20)) Output Gap Response to an Interest Rate Shock
med_gr_s_8.addtext(t, font("Times",20)) Inflation Response to an Interest Rate Shock
med_gr_s_9.addtext(t, font("Times",20)) Interest Rate Response to an Interest Rate Shock

delete(noerr) med_gr_s_all

graph med_gr_s_all.merge med_gr_s_1 med_gr_s_2 med_gr_s_3 med_gr_s_4 med_gr_s_5 med_gr_s_6 med_gr_s_7 med_gr_s_8 med_gr_s_9
med_gr_s_all.addtext(t, font("Times",30)) Cho and Morano (2005) Sign Restriction Graph (MT Response) Using SRC Approach

'Create the "sign restriction" graph

delete(noerr) gr_* gp_*

smpl @last - !n + 1 @last

graph gr_1.line {%vnames_1}
gr_1.options -legend
gr_1.addtext(t, font("Times",20)) Output Gap Response to a Demand Shock
gr_1.axis(l) zeroline 'range(-0.1, 1.0)
group gp_1.add {%vnames_1}

graph gr_2.line {%vnames_2}
gr_2.options -legend
gr_2.addtext(t, font("Times",20)) Inflation Response to a Demand Shock
gr_2.axis(l) zeroline 'range(-1.0, 0.1)
group gp_2.add {%vnames_2}

graph gr_3.line {%vnames_3}
gr_3.options -legend
gr_3.addtext(t, font("Times",20)) Interest Rate Response to a Demand Shock
gr_3.axis(l) zeroline 'range(-0.80, 0.1)
group gp_3.add {%vnames_3}

graph gr_4.line {%vnames_4}
gr_4.options -legend
gr_4.addtext(t, font("Times",20)) Output Gap Response to a Cost-Push Shock
gr_4.axis(l) zeroline 'range(-0.1, 1.0)
group gp_4.add {%vnames_4}

graph gr_5.line {%vnames_5}
gr_5.options -legend
gr_5.addtext(t, font("Times",20)) Inflation Response to a Cost-Push Shock
gr_5.axis(l) zeroline 'range(-0.1, 1.0)
group gp_5.add {%vnames_5}

graph gr_6.line {%vnames_6}
gr_6.options -legend
gr_6.addtext(t, font("Times",20)) Interest Rate Response to a Cost-Push Shock
gr_6.axis(l) zeroline 'range(-1.0, 0.1)
group gp_6.add {%vnames_6}

graph gr_7.line {%vnames_7}
gr_7.options -legend
gr_7.addtext(t, font("Times",20)) Output Gap Response to an Interest Rate Shock
gr_7.axis(l) zeroline 'range(-0.1, 1.0)
group gp_7.add {%vnames_7}

graph gr_8.line {%vnames_8}
gr_8.options -legend
gr_8.addtext(t, font("Times",20)) Inflation Response to an Interest Rate Shock
gr_8.axis(l) zeroline 'range(-0.1, 1.0)
group gp_8.add {%vnames_8}

graph gr_9.line {%vnames_9}
gr_9.options -legend
gr_9.addtext(t, font("Times",20)) Interest Rate Response to an Interest Rate Shock
gr_9.axis(l) zeroline 'range(-0.1, 1.0)
group gp_9.add {%vnames_9}

graph gr_all.merge gr_1 gr_2 gr_3 gr_4 gr_5 gr_6 gr_7 gr_8 gr_9
gr_all.addtext(t, font("Times",30)) Cho and Morano(2005) Sign Restriction Graph Using SRC

smpl @last - !n + 1 @last

show gr_all

!rr = (counter/!x_save)*100

@uiprompt("Retention rate = " + @str(!rr) + " (percent).")

@uiprompt("Press enter for the MT graph, SRC approach.") 

show med_gr_s_all


