'This program requires EViews 8 or higher to run...

subroutine local svar_impulses(scalar length, matrix shocks, matrix a0, matrix vlags, scalar neqs, scalar nlags, matrix result, scalar bool_accum)

	'Calculates the impulse response function of a structural VAR.

	'Input: length; scalar; number of inpulse responses to calculate, starting from period 0.

	'Input: a0; matrix; reflects the contemporaneous relationships in the system;it need not be triangular. 

	'Input: structural shocks; typically (though need not be) a diagonal matrix

	'One can scale the shocks (e.g. by the variance of the structural error) by placing the scalars on the diagonal of the "shock" matrix, which should be a square matrix, and have zeros on all the off diagonal elements.

	'Input: VLAGS; matrix; contains the lags in the structural VAR, arranged as lag1~lag2...lagN.  It can contain zero elements (representing zero restrictions on the SVAR lags). 

	'Rows of VLAGS = # of equations in the system; cols = #lag_length * rows of vlags.

	'Input: neqs; scalar; number of equations in the SVAR.

	'Input: nlags; scalar; maximum order of lags in the SVAR

	'Input: result; matrix; on completion, contains the impulse response functions.
	
	'Input: bool_accum  = 1 then return accumulated impulse response functions.

	matrix a0_inverse = @inverse(a0)

	for !j = 1 to nlags
		matrix phi{!j} = a0_inverse*@subextract(vlags, 1, 1+(!j-1)*neqs, neqs, neqs*!j)  'Note: these are lags of the VAR.
	next

	matrix(neqs, neqs) AAAA0 = @identity(neqs) 

	for !i=1 to nlags
		matrix(neqs, neqs) AAAA{!i} = 0
		for !j=1 to !i
			!imj = !i - !j
			AAAA{!i} = AAAA{!i} + (Phi{!j} * AAAA{!imj})
		next
	next

	for !i=nlags+1 to length-1
		matrix(neqs, neqs) AAAA{!i} = 0
		for !j=1 to nlags
			!imj = !i - !j
			AAAA{!i} = AAAA{!i} + (Phi{!j} * AAAA{!imj})
		next
	next

	matrix AAAA0Sigma = a0_inverse * shocks

	for !i=1 to length-1
		matrix  AAAA{!i}Sigma = AAAA{!i}*AAAA0Sigma
	next

	'Move the impulse responses to the "result" matrix, ordered by the response variables

	matrix(length, neqs*neqs) result = na

	for !hh=0 to length-1
		for !i = 1 to neqs
			for !j=1 to neqs
				result(!hh+1, ((!i-1)*neqs)+!j) =AAAA{!hh}Sigma(!i, !j)
			next
		next
	next

	'Accumulate the responses if bool_accum = 1

	if (bool_accum = 1) then
		for !hh = 2 to length
			for !j = 1 to (neqs^2)
				result(!hh,!j) = result(!hh-1,!j) + result(!hh,!j) 
			next
		next
	endif

endsub

subroutine local var_matrices(vector in_b, matrix ahat_svar, matrix alags, matrix f, matrix beta)

	'This subroutine is where the user sets up the matrices needed to estimate the SVAR Ay(t) = By(t-1)+Fx(t) + e(t) subject to arbitrary restrictions on B, F and the covariance matrix of e(t) 
	'B contains the coefficients of the lag endogenous variables. Each row of B represents a single equation.
	'The ordering is lag_1 (for all endogenous variables), followed by lag_2 (for all endogeneous variables), and so on...
	'X denotes the exogenous variable list. 
	'Each row of F contains coefficients on the exogenous variables. 
	'	
	'p_count (an internal scalar) keeps a running count of the parameter estimates that have been allocated/need to be estimated.
	'Variance matrix of  e(t) can be non-diagonal in general, but for this example we are assuming the matrix is diagonal (so we are estimating an SVAR).  
     'The user then needs to allocate these parameters to the A, B, F and covariance matrices. 
 
	'This allocation process can be automated if there no restrictions; more care needed if there are restrictions.

	!nvars = @rows(alags)  'Number of variables in the VAR
	!nlags = @columns(alags)/!nvars 'Number of lags in the VAR

	'Setup the contemporaneous A matrix, allowing for the long-run constraints.

	ahat_svar = 0

	!p = 0

	ahat_svar(1,1) = 1.0  'Diagonal terms
	ahat_svar(2,2) = 1.0
	ahat_svar(3,3) = 1.0
	ahat_svar(4,4) = 1.0
	ahat_svar(5,5) = 1.0

	ahat_svar(1,2) = in_b(!p+2) + in_b(!p+7) + in_b(!p+12) + in_b(!p+17)
	ahat_svar(1,3) = in_b(!p+3) + in_b(!p+8) + in_b(!p+13) + in_b(!p+18)
	ahat_svar(1,4) = in_b(!p+4) + in_b(!p+9) + in_b(!p+14) + in_b(!p+19)
	ahat_svar(1,5) = in_b(!p+5) + in_b(!p+10) + in_b(!p+15) + in_b(!p+20)

	'Now the lags. There are 100 lag coefficients to estimate: 1 to 100 

	!p_lags = !p+1

	for !i= 1 to 5
		for !j = 1 to 20
			alags(!i,!j) = in_b(!p_lags)
			!p_lags = !p_lags + 1
		next
	next

	'Count ends at 100 for the lagged coefficients+contemporaneous coeffcieints.  100 lagged coefficients plus 0 coefficients in the A matrix.

	'Now define the structure of the exogeneous coefficient matrix here; a constant term can be added to the VAR here (for example);  exogenous variables, etc, can be excluded from specific equations by using appropriate zero restrictions.

	'F will have a dimension of F(!nvars, !number_of_exogenous_variables).

	'Normally, the first column of F will represent the coefficients on the constant term. 

	f(1,1) = in_b(!p_lags) 'Constant terms;  !p_lags = 101
	f(2,1) = in_b(!p_lags+1)
	f(3,1) = in_b(!p_lags+2)
	f(4,1) = in_b(!p_lags+3)
	f(5,1) = in_b(!p_lags+4)

	'Now define the structure of the covariance matrix.   We will assume a diagonal covariance matrix.

	matrix beta = 0

	beta(1,1) = in_b(!p_lags+5) '  Variance of the first structural error

	beta(2,1) = in_b(!p_lags+6) 
	beta(2,2) = in_b(!p_lags+7)

	beta(3,1) = in_b(!p_lags+8)
	beta(3,2) = in_b(!p_lags+9)
	beta(3,3) = in_b(!p_lags+10)

	beta(4,1) = in_b(!p_lags+11)
	beta(4,2) = in_b(!p_lags+12)
	beta(4,3) = in_b(!p_lags+13)
	beta(4,4) = in_b(!p_lags+14)

	beta(5,1) = in_b(!p_lags+15)
	beta(5,2) = in_b(!p_lags+16)
	beta(5,3) = in_b(!p_lags+17)
	beta(5,4) = in_b(!p_lags+18)
	beta(5,5) = in_b(!p_lags+19)

	beta(1,2) = beta(2,1)
	beta(1,3) = beta(3,1)
	beta(1,4) = beta(4,1)
	beta(1,5) = beta(5,1)

	beta(2,3) = beta(3,2)
	beta(2,4) = beta(4,2)
	beta(2,5) = beta(5,2)

	beta(3,4) = beta(4,3)
	beta(3,5) = beta(5,3)

	beta(4,5) = beta(5,4)

'120 coefficients in total

endsub

subroutine loglikelihood(series logl, vector in_b, matrix y, matrix ylags, matrix x)

	'logmode all

' 	Estimates a structural VAR by maximizing the likelihood function using OPTIMIZE().  

'	Note: on entry, the vector in_b should contain the latest parameter estimates.

	!nvars = @columns(y)

	!n_ex  = @columns(x)

	!nlags = @columns(ylags)/!nvars

	matrix(!nvars, !nvars) x0 = 0

	matrix(!nvars,!nvars*!nlags) vlags = 0

	matrix(!nvars, !n_ex) f = 0

	matrix(!nvars,!nvars) beta = 0

	matrix(!nvars,!nvars) ahat_svar = 0

	call var_matrices(in_b, ahat_svar, vlags, f, beta)

	'Now construct the residual matrix for the current parameter estimates...and then calculate the value of the likelihood function.

	matrix r = (y*@transpose(ahat_svar)-(ylags*@transpose(vlags)))-(x*@transpose(f)) 'VAR residuals...notice A (or x0) is assumed to be the identity matrix.

	!det_x0 = @det(ahat_svar) ' = 1 if A matrix is triangular.

	matrix sigma = beta'*@transpose(beta) 'It's actually beta*@transpose(beta), but this is OK, since we then do not need to unravel the coefficients.

	!det_sigma = @det(sigma)

	matrix i_sigma = @inverse(sigma)

	vector(@rows(r),1) v

	!t = @rows(r)

	for !j = 1 to !t
		rowvector y_t = @rowextract(r,!j)
		!pp = (y_t*i_sigma*@transpose(y_t))
		v(!j,1) = {!pp}
	next

	mtos(v,vv)  'Need to convert the vector to a series...this guarantees that we get standard errors in the end

	'Likelihood function; note that the new optimizer allows the elements (e.g., determinant of A) to be a function of the unknown parameters...

	'Note: pi() = @acos(-1)

	!pi = @acos(-1)

	logl = -(!nvars/2.0)*log(2*!pi) - (1/2)*log(!det_sigma) + log(!det_x0) - (0.5*vv)

endsub

'program uses the data in galidusa.wf1

'This code replicates Gali (1999): "Technology, Employment and the Business Cycle: Do Technology Shocks Explain Aggregate Fluctuations
' AER (1999), pages 249-271.

pageselect Data

'logmode all

!m = 5				'Number of equations
!var_lags  = 4		'Number of lags
!nx = 28				'Number of periods required in the impulse vector
!b_iters = 1     		'Number of times to change the underlying VAR parameter estimates during the bootstrap; (not used for this example).
!r_iters = 1    		'Number of iterations for a residual bootstrap, given the current estimates for the SVAR; (not used for this example)
!bm_method = 1	'Set = 1 to bootstrap the residuals from the N(0,1) distribution
!debug = -1			'False = < 0; True >= 0; if false, suppresses most error messages.
!npars = 120	     	'Total number of parameters to estimate in the system; (100 lags; 5 constants; and 15 in the covariance matrix)
series l1				'Will eventually contain the value of the likelihood at each observation.

'Create the data matrices

smpl @all

delete(noerr) vdata* xm*

'smpl 1947q1 1994q4

smpl 1960q2 1994q4

group vdata.add dprodh dhours dinf ec1 ec2

group xdata.add c

stom(xdata,xm) 

stom(vdata,ym)

delete(noerr) ylags*

group ylags

ylags.add dprodh(-1) dhours(-1) dinf(-1) ec1(-1) ec2(-1) dprodh(-2) dhours(-2) dinf(-2) ec1(-2) ec2(-2) dprodh(-3) dhours(-3) dinf(-3) ec1(-3) ec2(-3) dprodh(-4) dhours(-4) dinf(-4) ec1(-4) ec2(-4) 

'Notice the ordering of the lags: all of -1 first; followed by all of -2 second, and so on.

stom(ylags, ym_lags)

smpl 1960q2 1994q4

var galitech.ls 1 4  dprodh dhours dinf ec1 ec2  @ c

show galitech

vector(!npars) fe = 0

matrix sigma= @identity(!m) 'This will eventually contain the variances of the structural errors

matrix ui = galitech.@coefmat

!z =0

fe(1,1) = ui(1,1)
fe(2,1) = ui(1,2)
fe(3,1) = ui(1,3)
fe(4,1) = ui(1,4)
fe(5,1) = ui(1,5)

fe(6,1) = ui(5,1)
fe(7,1) = ui(5,2)
fe(8,1) = ui(5,3)
fe(9,1) = ui(5,4)
fe(10,1) = ui(5,5)

fe(11,1) = ui(9,1)
fe(12,1) = ui(9,2)
fe(13,1) = ui(9,3)
fe(14,1) = ui(9,4)
fe(15,1) = ui(9,5)
'
fe(16,1) = ui(13,1)
fe(17,1) = ui(13,2)
fe(18,1) = ui(13,3)
fe(19,1) = ui(13,4)
fe(20,1) = ui(13,5)

fe(21,1) = ui(17,1)
fe(22,1) = ui(17,2)
fe(23,1) = ui(17,3)
fe(24,1) = ui(17,4)
fe(25,1) = ui(17,5)

fe(26,1) = ui(2,1)
fe(27,1) = ui(2,2)
fe(28,1) = ui(2,3)
fe(29,1) = ui(2,4)
fe(30,1) = ui(2,5)

fe(31,1) = ui(6,1)
fe(32,1) = ui(6,2)
fe(33,1) = ui(6,3)
fe(34,1) = ui(6,4)
fe(35,1) = ui(6,5)

fe(36,1) = ui(10,1)
fe(37,1) = ui(10,2)
fe(38,1) = ui(10,3)
fe(39,1) = ui(10,4)
fe(40,1) = ui(10,5)
'
fe(41,1) = ui(14,1)
fe(42,1) = ui(14,2)
fe(43,1) = ui(14,3)
fe(44,1) = ui(14,4)
fe(45,1) = ui(14,5)

fe(46,1) = ui(18,1)
fe(47,1) = ui(18,2)
fe(48,1) = ui(18,3)
fe(49,1) = ui(18,4)
fe(50,1) = ui(18,5)

fe(51,1) = ui(3,1)
fe(52,1) = ui(3,2)
fe(53,1) = ui(3,3)
fe(54,1) = ui(3,4)
fe(55,1) = ui(3,5)

fe(56,1) = ui(7,1)
fe(57,1) = ui(7,2)
fe(58,1) = ui(7,3)
fe(59,1) = ui(7,4)
fe(60,1) = ui(7,5)

fe(61,1) = ui(11,1)
fe(62,1) = ui(11,2)
fe(63,1) = ui(11,3)
fe(64,1) = ui(11,4)
fe(65,1) = ui(11,5)
'
fe(66,1) = ui(15,1)
fe(67,1) = ui(15,2)
fe(68,1) = ui(15,3)
fe(69,1) = ui(15,4)
fe(70,1) = ui(15,5)

fe(71,1) = ui(19,1)
fe(72,1) = ui(19,2)
fe(73,1) = ui(19,3)
fe(74,1) = ui(19,4)
fe(75,1) = ui(19,5)

fe(76,1) = ui(4,1)
fe(77,1) = ui(4,2)
fe(78,1) = ui(4,3)
fe(79,1) = ui(4,4)
fe(80,1) = ui(4,5)

fe(81,1) = ui(8,1)
fe(82,1) = ui(8,2)
fe(83,1) = ui(8,3)
fe(84,1) = ui(8,4)
fe(85,1) = ui(8,5)

fe(86,1) = ui(12,1)
fe(87,1) = ui(12,2)
fe(88,1) = ui(12,3)
fe(89,1) = ui(12,4)
fe(90,1) = ui(12,5)
'
fe(91,1) = ui(16,1)
fe(92,1) = ui(16,2)
fe(93,1) = ui(16,3)
fe(94,1) = ui(16,4)
fe(95,1) = ui(16,5)

fe(96,1) = ui(20,1)
fe(97,1) = ui(20,2)
fe(98,1) = ui(20,3)
fe(99,1) = ui(20,4)
fe(100,1) = ui(20,5)

fe(101,1) = ui(21,1) 'Constant terms
fe(102,1) = ui(21,2)
fe(103,1) = ui(21,3)
fe(104,1) = ui(21,4)
fe(105,1) = ui(21,5)

matrix temp = galitech.@residcov

fe(!npars-14,1) = temp(1,1) 

fe(!npars-13,1) = temp(2,1) 
fe(!npars-12,1) = temp(2,2) 

fe(!npars-11,1) = temp(3,1)
fe(!npars-10,1) = temp(3,2)
fe(!npars-9,1) = temp(3,3)

fe(!npars-8,1) = temp(4,1)
fe(!npars-7,1) = temp(4,2)
fe(!npars-6,1) = temp(4,3)
fe(!npars-5,1) = temp(4,4)

fe(!npars-4,1) = temp(5,1)
fe(!npars-3,1) = temp(5,2)
fe(!npars-2,1) = temp(5,3)
fe(!npars-1,1) = temp(5,4)
fe(!npars,1) = temp(5,5)

matrix fe_starting_values = fe

show fe

matrix rcov = galitech.@residcov

optimize(max=1,finalh=mlhess,hess=bfgs,deriv=auto,m=100000000,c=1e-10) loglikelihood(l1,fe,ym,ym_lags,xm)

vector ml_se = @sqrt(@getmaindiagonal(-@inverse(mlhess)))  'ML estimates of the standard errors; closest to .'observed hessian option in EViews.

scalar llike_value = @sum(l1)  'Value of the likelihood function.

delete(noerr) rr

matrix(!npars,3) rr

colplace(rr,fe,1)
colplace(rr,ml_se,2)
matrix(!npars,1) vt = @ediv(fe,ml_se)
colplace(rr,vt,3)

delete vt


matrix a0 = @identity(!m)

matrix(!m,!m*!var_lags) a1 = 0

'a1 will contain the lag parameter estimates.  

matrix(!m,2) c1 = 0 ' Exogenous variable coefficients

matrix shocks = @identity(!m)

call var_matrices(fe, a0, a1, c1,sigma)

matrix a0_inv = @inverse(a0)

'Shock matrix allows the user to scale the initial shock to the VAR.  The impulse function routine uses a default shock of 1.0

matrix shocks = @identity(!m)

shocks(1,1) = sigma(1,1)^0.5
shocks(2,2) = sigma(2,2)^0.5
shocks(3,3) = sigma(3,3)^0.5
shocks(4,4) = sigma(4,4)^0.5

matrix(!nx,!m) impulses = na

'Now calculate the impulse response function for the first parameter estimates (initial convergence..

call svar_impulses(!nx, shocks, a0, a1, !m, !var_lags, impulses,1)

delete(noerr) impulse_acc*

matrix impulse_acc = impulses

show impulse_acc

'matrix diff = ahat_svar-fiml_ahat

'show diff


