function mGa = ComputeMeanLag_svec_var(mGamma)
% /*
% **	Purpose: Compute mean lag matrix needed for computation of long
% **			 long run impact matrix C(1).
% **
% **	Imput:   mGamma:  K x K(p-1) matrix coefficient matrix of lagged differences 
% **			    	  \Gamma = \Gamma_1 ... \Gamma_{p-1}  
% **
% **	Output:  mGa:    I_K - Gamma_1 - ... - \Gamma_{p-1}
% **
% **                  equals I_K if no lagged differences are included 
% */

K 	= size(mGamma,1);
p 	= size(mGamma,2)/K;
mGa = eye(K);

i = 1;
while i <= p
    mGa = mGa - mGamma(:,(i-1)*K+1:i*K);
    i = i +1;
end
