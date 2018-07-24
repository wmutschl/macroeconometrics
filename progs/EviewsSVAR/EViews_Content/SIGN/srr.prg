smpl @all

pageselect chomoreno

delete(noerr) impulse_set_* counter* gp_* gr_* i*_* svar_impulse* sign_mat* inst* e1 e2 e3 a_* vv*

!nvars = 3         	'Number of variables in the base VAR

!n = 32 			     'Length of a single impulse response function

!window = !n		'Window: criterion: first eight values of the impulse response vector must have the correct sign

!j   = (!nvars^2) 	'Number of impulse response vectors in a single run

!niters = 10000 'Total number of impulse response function sets to generate; number taken will typically be much smaller.

rndseed 123456789

for !s = 1 to !j
	matrix(!n,!niters) impulse_set_{!s}=na  'This will ultimately store all the valid impuse response vectors
	%vnames_{!s} = "" 'Keeps track of the names of these vectors...
next

delete(noerr) sign_matrix ch

var ch.ls 1 2 gap infl ff 'Create the VAR

%smpl = ch.@smpl

scalar counter = 0  'This will keep track of the valid impulse count for each shock.

matrix(!nvars, !nvars) q12 =  0 

matrix(!nvars, !nvars) q13 = 0

matrix(!nvars, !nvars) q23 = 0

matrix(!nvars, !nvars) pattern_a = 0

matrix(!nvars, !nvars) pattern_b = 0

scalar success = 0

'Gap equation...

equation ls_gap.ls gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c 

pattern_b(1,1) = ls_gap.@se

pattern_a(1,1) = 1.0

'Inflation equation...

equation ls_infl.ls infl gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2)  c 

pattern_a(2,1) = -ls_infl.@coefs(1)

pattern_a(2,2) = 1.0

pattern_a(2,3) = 0.0

pattern_b(2,2) = ls_infl.@se

'Interest rate equation

equation ls_ff.ls ff infl gap gap(-1) gap(-2) infl(-1) infl(-2) ff(-1) ff(-2) c 

pattern_a(3,2) = -ls_ff.@coefs(1)

pattern_a(3,1) = -ls_ff.@coefs(2)

pattern_a(3,3) = 1

pattern_b(3,3) = ls_ff.@se

matrix inv_pattern_b = @inverse(pattern_b)

!pi = @acos(-1)

for !x = 1 to !niters 'Iteration loop

	'Generate a new set of identifying parameters...

	matrix q12 =  0 

	matrix q13 = 0

	matrix q23 = 0

	!p1 = @runif(0,1)*!pi

	!p2 = @runif(0,1)*!pi

	!p3 = @runif(0,1)*!pi

	q12(1,1) = @cos(!p1)

	q12(1,2) = -@sin(!p1)

	q12(2,1) = @sin(!p1)

	q12(2,2) = @cos(!p1)

	q12(3,3) = 1

	q13(1,1) = @cos(!p2)

	q13(1,3) = -@sin(!p2)

	q13(3,1) = @sin(!p2)

	q13(3,3) = @cos(!p2)

	q13(2,2) = 1

	q23(2,2) = @cos(!p3)
	
	q23(2,3) = -@sin(!p3)

	q23(3,2) = @sin(!p3)

	q23(3,3) = @cos(!p3)

	q23(1,1) = 1

	matrix q = q12*q13*q23

	matrix a_inv = @inverse(q*inv_pattern_b*pattern_a)

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

	'Demand

	for !j = 1 to 3

		!test_p = ((a_inv(1,!j) >= 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) >= 0))

		if !test_p > 0 then 

			!pd_count = !pd_count + 1

			!d_column = !j

		endif

		!test_n = ((a_inv(1,!j) < 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) < 0))

		if !test_n > 0 then 

			!nd_count = !nd_count + 1

			!d_column = -!j

		endif

	next

	if ((!pd_count+!nd_count) = 1) then 

		!t1 = 1

	endif

	'Inflation

	for !j = 1 to 3

		!test_p = ((a_inv(1,!j) < 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) >= 0))

		if !test_p > 0 then 

			!pp_count = !pp_count + 1

			!p_column = !j

		endif

		!test_n = ((a_inv(1,!j) >= 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) < 0))

		if !test_n > 0 then 

			!np_count = !np_count + 1

			!p_column = -!j

		endif

	next

	if ((!pp_count+!np_count) = 1) then 

		!t2 = 1

	endif

	'Interest rate shock...

	for !j = 1 to 3

		!test_p = (a_inv(1,!j) < 0) and (a_inv(2,!j) < 0) and (a_inv(3,!j) >= 0)

		if !test_p > 0 then 

			!pr_count = !pr_count + 1
	
			!r_column = !j

		endif

		!test_n = (a_inv(1,!j) >= 0) and (a_inv(2,!j) >= 0) and (a_inv(3,!j) < 0)

		if !test_n > 0 then 

			!nr_count = !nr_count + 1

			!r_column = -!j

		endif

	next

	if ((!pr_count+!nr_count) = 1) then 

		!t3 = 1

	endif

	if ((!t1+!t2 +!t3)>2.5) then
	
		success = success + 1

'		Re-order the columns so that demand shock somes first, cost-push is second, and interest rate is third.

		matrix final_a = a_inv

		vector v1 = @subextract(a_inv,1,@abs(!d_column),!nvars,@abs(!d_column))

		if !d_column < 0 then
		 	v1 = -v1 'Negative demand shock; make it positive and store it in column 1
		endif

		colplace(final_a,v1,1)

		vector v2 = @subextract(a_inv,1,@abs(!p_column),!nvars,@abs(!p_column))

		if !p_column < 0 then
			v2 = -v2  'Negative cost push shock; make it positive and store it in column 2
		endif

		colplace(final_a,v2,2)

		vector v3 = @subextract(a_inv,1,@abs(!r_column),!nvars,@abs(!r_column))

		if !r_column < 0 then
			v3 = -v3 'Negative interest rate shock; make it positive and store it in column 3
		endif

		colplace(final_a,v3,3)

		delete(noerr) ch_svar_imp* 'Remove the previous impulse response function set

'		Calculate the new impulse response function...assuming the variance of the errors are unitary (across all errors)

'		Impulse reponse functions are ordered as gap w.r.t gap, infl w.r.t gap, ff with respect to gap
'																gap w.r.t infl,   infl w.r.t. infl; ff w.r.t. infl
'																gap w.r.t ff, infl w.r.t ff, ff w.r.t ff	
'							
		freeze(mode=overwrite,impgraph) ch.impulse(!n,smat=ch_svar_imp,se=a,imp=user,fname=final_a) gap infl ff @ gap infl ff @ gap infl ff
		
		counter = counter+1  'Update the valid counter index 
	
		%tt = @str(counter)

		%tt ="Success count = " + @str(counter) + " models; " + @str((counter/!niters)*100) + " percent." 

		statusline {%tt}

		for !v = 1 to !nvars^2

			delete(noerr) p1

			smpl @last - !n + 1 @last

			!t = counter

			vector p1 = @columnextract(ch_svar_imp,!v)

			colplace(impulse_set_{!v},p1,!t)

			series g{!v}_{!t} = 0 'In addition to the impulse_set store, create a series for the valid impulse response function

			mtos(p1,g{!v}_{!t}) 'Save the valid impulse to a series (for graphing purposes)
                                                                       		
			%vnames_{!v} = %vnames_{!v} + " " + "g"+@str({!v})+"_"+@str({!t}) 'Keep track of the valid series names (each is an impulse response)...

		next

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

!oo = 10^14

for !cc = 1 to counter

	matrix r_pk = @rowextract(pk,!cc)

	!oox = r_pk*@transpose(r_pk)

	xmin(!cc,1) = !oox

	if xmin(!cc,1) <= !oo then 

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
	delete(noerr) median_imp_g_{!v}
	series median_imp_g_{!v} = 0
	mtos(c_temp, median_imp_g_{!v})
next

'Create the "sign restriction" graph

delete(noerr) gr_* gp_*

smpl @last - !n + 1 @last

for !v = 1 to !nvars^2
	delete(noerr) med_gr_{!v}
	graph med_gr_{!v}.line median_imp_g_{!v}
	med_gr_{!v}.options -legend
	med_gr_{!v}.axis(l) zeroline
next

med_gr_1.addtext(t, font("Times",20)) Output Gap Response to a Demand Shock
med_gr_2.addtext(t, font("Times",20)) Inflation Response to a Demand Shock
med_gr_3.addtext(t, font("Times",20)) Interest Rate Response to a Demand Shock
med_gr_4.addtext(t, font("Times",20)) Output Gap Response to a Cost-Push Shock
med_gr_5.addtext(t, font("Times",20)) Inflation Response to a Cost-Push Shock
med_gr_6.addtext(t, font("Times",20)) Interest Rate Response to a Cost-Push Shock
med_gr_7.addtext(t, font("Times",20)) Output Gap Response to an Interest Rate Shock
med_gr_8.addtext(t, font("Times",20)) Inflation Response to an Interest Rate Shock
med_gr_9.addtext(t, font("Times",20)) Interest Rate Response to an Interest Rate Shock

delete(noerr) med_gr_all

graph med_gr_all.merge med_gr_1 med_gr_2 med_gr_3 med_gr_4 med_gr_5 med_gr_6 med_gr_7 med_gr_8 med_gr_9
med_gr_all.addtext(t, font("Times",30)) Cho and Morano (2005) Sign Restriction Graph (Median Response) Using Givens Approach

graph gr_g1.line {%vnames_1}
gr_g1.options -legend
gr_g1.addtext(t, font("Times",20)) Output Gap Response to a Demand Shock
gr_g1.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g1.add {%vnames_1}

graph gr_g2.line {%vnames_2}
gr_g2.options -legend
gr_g2.addtext(t, font("Times",20)) Inflation Response to a Demand Shock
gr_g2.axis(l) zeroline 'range(-1.0, 0.1)
group gp_g2.add {%vnames_2}

graph gr_g3.line {%vnames_3}
gr_g3.options -legend
gr_g3.addtext(t, font("Times",20)) Interest Rate Response to a Demand Shock
gr_g3.axis(l) zeroline 'range(-0.80, 0.1)
group gp_g3.add {%vnames_3}

graph gr_g4.line {%vnames_4}
gr_g4.options -legend
gr_g4.addtext(t, font("Times",20)) Output Gap Response to a Cost-Push Shock
gr_g4.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g4.add {%vnames_4}

graph gr_g5.line {%vnames_5}
gr_g5.options -legend
gr_g5.addtext(t, font("Times",20)) Inflation Response to a Cost-Push Shock
gr_g5.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g5.add {%vnames_5}

graph gr_g6.line {%vnames_6}
gr_g6.options -legend
gr_g6.addtext(t, font("Times",20)) Interest Rate Response to a Cost-Push Shock
gr_g6.axis(l) zeroline 'range(-1.0, 0.1)
group gp_g6.add {%vnames_6}

graph gr_g7.line {%vnames_7}
gr_g7.options -legend
gr_g7.addtext(t, font("Times",20)) Output Gap Response to an Interest Rate Shock
gr_g7.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g7.add {%vnames_7}

graph gr_g8.line {%vnames_8}
gr_g8.options -legend
gr_g8.addtext(t, font("Times",20)) Inflation Response to an Interest Rate Shock
gr_g8.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g8.add {%vnames_8}

graph gr_g9.line {%vnames_9}
gr_g9.options -legend
gr_g9.addtext(t, font("Times",20)) Interest Rate Response to an Interest Rate Shock
gr_g9.axis(l) zeroline 'range(-0.1, 1.0)
group gp_g9.add {%vnames_9}

graph gr_all_g.merge gr_g1 gr_g2 gr_g3 gr_g4 gr_g5 gr_g6 gr_g7 gr_g8 gr_g9
gr_all_g.addtext(t, font("Times",30)) Cho and Morano (2005) Sign Restriction Graph Using Givens Approach

smpl @last - !n + 1 @last

show gr_all_g

!rr = (counter/!niters)*100

@uiprompt("Retention rate = " + @str(!rr) + " (percent), i.e., " + @str(counter) + " models.")

@uiprompt("Press enter for the median response graph, Givens approach.") 

show med_gr_all

