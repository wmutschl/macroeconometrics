function mC1 = ComputeC1_svec_var(alpha, beta, mGamma)
% /*
% **	Usage: mC1 = ComputeC1(alpha, beta, mGamma);
% **
% **	Purpose: Compute long run impact matrix for cointegrated VAR models
% **
% **	Input:		alpha:	 K x r matrix of loading coefficients
% **
% **				beta: 	 K x r matrix of cointegration coefficients related to stochastic variables
% **
% **				mGamma:  K x K(p-1) matrix coefficient matrix of lagged differences 
% **				  		 \Gamma = \Gamma_1 ... \Gamma_{p-1}  
% **
% **	Ouput: 		mC1:	 K x K matrix, reduced form long run impact matrix
% */
K = size(alpha,1);
if isempty(mGamma)			%/* check whether lags are present*/ 
    mGa = eye(K);
else
    mGa = ComputeMeanLag_svec_var(mGamma);
end
		
beta_o = null(beta');	%/* orthogonal complement of beta */
alpha_o = null(alpha');	%/* orthogonal complemt of alpha  */

mC1 = beta_o*inv(alpha_o'*mGa*beta_o)*alpha_o'; %/* total impact matrix C(1) */