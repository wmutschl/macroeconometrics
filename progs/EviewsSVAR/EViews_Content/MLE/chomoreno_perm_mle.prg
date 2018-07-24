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

	ahat_svar(1,1) = 1.0  'Diagonal terms
	ahat_svar(1,2) = in_b(5)+in_b(8)
	ahat_svar(1,3) = in_b(6)+in_b(9)

	ahat_svar(2,1) = in_b(1)
	ahat_svar(2,2) = 1.0
	ahat_svar(2,3) = in_b(2)

	ahat_svar(3,1) = in_b(3)
	ahat_svar(3,2) = ahat_svar(1,2)/ahat_svar(1,3)
	ahat_svar(3,3) = 1.0
	
	'show ahat_svar

	'Now the lags. There are 18 lag coefficients to estimate: 4 (inclusive) to 21

	!p_lags = 4

	for !i= 1 to 3
		for !j = 1 to 6
			alags(!i,!j) = in_b(!p_lags)
			!p_lags = !p_lags + 1
		next
	next

	'Count ends at 51 for the lagged coefficients+contemporaneous coeffcieints.  18 lagged coefficients plus 3 coefficients in the A0 matrix.

	'Now define the structure of the exogeneous coefficient matrix here; a constant term can be added to the VAR here (for example);  exogenous variables, etc, can be excluded from specific equations by using appropriate zero restrictions.

	'F will have a dimension of F(!nvars, !number_of_exogenous_variables).

	'Normally, the first column of F will represent the coefficients on the constant term. 

	f(1,1) = in_b(22) 'Constant terms
	f(2,1) = in_b(23)
	f(3,1) = in_b(24)

	'Now define the structure of the covariance matrix.   We will assume a diagonal covariance matrix.

	matrix beta = 0

	beta(1,1) = in_b(25) '  Variance of the first structural error
	beta(2,2) = in_b(26) '  ...
	beta(3,3) = in_b(27)

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

pageselect Chomoreno

'logmode all

!m = 3				'Number of equations
!var_lags  = 2		'Number of lags
!nx = 28				'Number of periods required in the impulse vector
!b_iters = 1     		'Number of times to change the underlying VAR parameter estimates during the bootstrap; (not used for this example).
!r_iters = 1    		'Number of iterations for a residual bootstrap, given the current estimates for the SVAR; (not used for this example)
!bm_method = 1	'Set = 1 to bootstrap the residuals from the N(0,1) distribution
!debug = -1			'False = < 0; True >= 0; if false, suppresses most error messages.
!npars = 27      	'Total number of parameters to estimate in the system; (18 lags; 3 constants; 3 contemporaneous; and 3 in the covariance matrix)
series l1				'Will eventually contain the value of the likelihood at each observation.

'Create the data matrices

smpl @all

delete(noerr) vdata* xm*

smpl 1981q1 2001q2

group vdata.add dgap ff infl

group xdata.add c

smpl 1981q4 2000q1

stom(xdata,xm)

stom(vdata,ym)

delete(noerr) ylags*

group ylags

ylags.add dgap(-1) ff(-1) infl(-1) dgap(-2) ff(-2) infl(-2) 

'Notice the ordering of the lags: all of -1 first; followed by all of -2 second, and so on.

stom(ylags, ym_lags)

var chomorperm.ls 1 2 dgap ff infl @ c

show chomorperm

'The constrained system imposes summing up constraints and a cross-equation restriction ....

vector(!npars) fe = 0.1

matrix sigma=@identity(!m) 'This will eventually contain the variances of the structural errors

chomoreno_sys.fiml(conv=1e-10,rcov=diag)

matrix ahat_fiml = @identity(!m)

ahat_fiml(1,3) = (c(3)+c(4))
ahat_fiml(1,2) = (c(5)+c(6))
ahat_fiml(2,1) = -c(16)
ahat_fiml(2,3) = -c(17)
ahat_fiml(3,1) = -c(25)
ahat_fiml(3,2) = (c(5)+c(6))/(c(3)+c(4))

show ahat_fiml

'ahat_svar = a_start
'
'show ahat_svar
'
'fe(1,1) = ahat_svar(2,1)  'a_start(2,1) 'oil
'
'fe(2,1) = ahat_svar(3,1) '0 'a_start(3,1) '0'-c(62) 'inflation
'
'fe(3,1) = ahat_svar(3,2) 'a_start(3,2) '-c(63)  'Oil
'fe(4,1) = ahat_svar(3,4) '0' a_start(3,4) '0'-c(64)  'gdp
'fe(5,1) = ahat_svar(4,1) 
'fe(6,1) = ahat_svar(4,2)
'
'show fe
'
'fe(7,1) = c(1) 'Lag 1 (First Equation)
'fe(8,1) = c(4)
'fe(9,1) = c(7)
'fe(10,1) = c(10)
'
'fe(11,1) = c(2) 'Lag 2  (First Equation)
'fe(12,1) = c(5)
'fe(13,1) = c(8)
'fe(14,1) = c(11)
'
'fe(15,1) = c(3) 'Lag 3  (First Equation)
'fe(16,1) = c(6)
'fe(17,1) = c(9)
'fe(18,1) = c(12)
'
'fe(19,1) = c(15) 'Lag 1 (Second Equation)
'fe(20,1) = c(18)
'fe(21,1) = c(21)
'fe(22,1) = c(24)
'
'fe(23,1) = c(16) 'Lag 2
'fe(24,1) = c(19)
'fe(25,1) = c(22)
'fe(26,1) = c(25)
'
'fe(27,1) = c(17) 'Lag 3
'fe(28,1) = c(20)
'fe(29,1) = c(23)
'fe(30,1) = c(26)
'
'fe(31,1) = c(29)  'Third Equation; lag 1
'fe(32,1) = c(32)
'fe(33,1) = c(35)
'fe(34,1) = c(38)
'
'fe(35,1) = c(30) 'Third equation, Lag 2
'fe(36,1) = c(33)
'fe(37,1) = c(36)
'fe(38,1) = c(39)
'
'fe(39,1) = c(31) 'Third equation, Lag 3
'fe(40,1) = c(34)
'fe(41,1) = c(37)
'fe(42,1) = c(40)
'
'fe(43,1) = c(43) 'Fourth Equation; lag 1
'fe(44,1) = c(46)
'fe(45,1) = c(49)
'fe(46,1) = c(52)
'
'fe(47,1) = c(44) 'Fourth Equation; lag 2
'fe(48,1) = c(47)
'fe(49,1) = c(50)
'fe(50,1) = c(53)
'
'fe(51,1) = c(45) 'Fourth Equation; lag 3
'fe(52,1) = c(48)
'fe(53,1) = c(51)
'fe(54,1) = c(54)
'
'fe(55,1) = c(13) 'Constant term
'fe(56,1) = c(27)
'fe(57,1) = c(41)
'fe(58,1) = c(55)
'
'fe(59,1) = c(14) 'Trend term
'fe(60,1) = c(28)
'fe(61,1) = c(42)
'fe(62,1) = c(56)

fe(25,1) = chomorperm.@residcov(1,1)
fe(26,1) = chomorperm.@residcov(2,2)
fe(27,1) = chomorperm.@residcov(3,3)

matrix fe_starting_values = fe

show fe

smpl 1981q4 2000q1

optimize(max=1,finalh=mlhess,hess=bfgs,deriv=high,m=10000000,c=1e-10) loglikelihood(l1,fe,ym,ym_lags,xm)

vector ml_se = @sqrt(@getmaindiagonal(-@inverse(mlhess)))  'ML estimates of the standard errors.

scalar llike_value = @sum(l1)  'Value of the likelihood function.

delete(noerr) rr

matrix(!npars,2) rr

colplace(rr,fe,1)
colplace(rr,ml_se,2)

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

matrix(!nx,!m) impulses = na

'Now calculate the impulse response function for the first parameter estimates (initial convergence..

call svar_impulses(!nx, shocks, a0, a1, !m, !var_lags, impulses,0)

delete(noerr) impulse_*

'matrix impulse_acc = impulses

'show impulse_acc

show impulses


