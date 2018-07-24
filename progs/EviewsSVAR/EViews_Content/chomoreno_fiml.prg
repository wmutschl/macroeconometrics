'requires chomoreno.wf1
pageselect Chomoreno
smpl @all
'Estimate the summative (reduced form) VAR
var chomoreno.ls 1 2 gap infl ff
chomoreno.results

'Estimate the same model using FIML, a lower-triangular identification
'scheme and a diagonal covariance matrix.

chomor_sys.fiml(covinfo=hessian, rcov=diag) 

'Build the contemporaneous A matrix

matrix ahat = @identity(3)
ahat(2,1) = -c(22)
ahat(3,1) = -c(23)
ahat(3,2) = -c(24)

'And the B matrix

matrix bhat = chomor_sys.@estcov
bhat(1,1) = bhat(1,1)^0.5
bhat(2,2) = bhat(2,2)^0.5
bhat(3,3) = bhat(3,3)^0.5

'Since the model is exactly identified, we can use the 
'summative model to calculate the impulse response
'functions.

matrix shocks = @inverse(ahat)*bhat
chomoreno.impulse(10,m,imp=user,se=a,fname=shocks)

