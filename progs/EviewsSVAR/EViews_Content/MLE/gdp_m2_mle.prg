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
	ahat_svar(2,2) = 1.0

	ahat_svar(1,2) = in_b(3) 'This forces c(1,2) = 0
	ahat_svar(2,1) = in_b(1)

	'Now the lags. There are 4 lags to estimate.

	!p_lags = 2

	for !i= 1 to !nvars
		for !j = 1 to (!nlags*!nvars)
			alags(!i,!j) = in_b(!p_lags)
			!p_lags = !p_lags + 1
		next
	next

	'Count ends at 5 for the lagged coefficients+contemporaneous coeffcieints.  4 lagged coefficients plus 1 coefficient in the @inverse(A0) matrix.

	'Now define the structure of the exogeneous coefficient matrix here; a constant term can be added to the VAR here (for example);  exogenous variables, etc, can be excluded from specific equations by using appropriate zero restrictions.

	'F will have a dimension of F(!nvars, !number_of_exogenous_variables).

	'Normally, the first column of F will represent the coefficients on the constant term. 

	f(1,1) = in_b(6) 'Constant terms
	f(2,1) = in_b(7)

	'Now define the structure of the covariance matrix.   We will assume a diagonal covariance matrix.

	matrix(2,2) beta = 0

	beta(1,1) = in_b(8) '  Variance of the first structural error
	beta(2,2) = in_b(9) '  ...

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

pageselect gdp_m2

'logmode all

!m = 2				'Number of equations
!var_lags  = 1		'Number of lags
!nx = 28				'Number of periods required in the impulse vector
!b_iters = 1     		'Number of times to change the underlying VAR parameter estimates during the bootstrap; (not used for this example).
!r_iters = 1    		'Number of iterations for a residual bootstrap, given the current estimates for the SVAR; (not used for this example)
!bm_method = 1	'Set = 1 to bootstrap the residuals from the N(0,1) distribution
!debug = -1			'False = < 0; True >= 0; if false, suppresses most error messages.
!npars = 9      		'Total number of parameters to estimate in the system; (4 lags; 2 constants; 1 contemporaneous; 2 in the covariance matrix)
series l1				'Will eventually contain the value of the likelihood at each observation.

'Create the data matrices

smpl @all

delete(noerr) vdata* xm* beta*

smpl 1981q3 2000q1

group vdata.add s_dgdp s_dm2

group xdata.add c

stom(xdata,xm)

stom(vdata,ym)

delete(noerr) ylags*

group ylags

ylags.add s_dgdp(-1) s_dm2(-1)

'Notice the ordering of the lags: all of -1 first; followed by all of -2 second, and so on.

stom(ylags, ym_lags)

var gdpm2.ls 1 1 s_dgdp s_dm2 @ c

show gdpm2

'The constrained system imposes summing up constraints and a cross-equation restriction ....

vector(!npars) fe =0

fe(!npars-1,1) = gdpm2.@residcov(1,1)
fe(!npars,1) = gdpm2.@residcov(2,2)

matrix sigma=@identity(!m) 'This will eventually contain the variances of the structural errors

matrix fe_starting_values = fe

show fe

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

matrix(!nx,!m) impulses = na

'Now calculate the impulse response function for the first parameter estimates (initial convergence..

call svar_impulses(!nx, shocks, a0, a1, !m, !var_lags, impulses,1)

delete(noerr) impulse_*

'matrix impulse_acc = impulses

'show impulse_acc

show impulses


