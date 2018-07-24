function [vgn,i_irt,i, noConvergence] = DecompCorrelation_svar_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model,imaxRetries)
% /*
% **	Purpose: compute structural decomposition of correlation matrix to find optimal start values
% **	
% **	Usage: {vgn,i_irt,i} = DecompCorrelation(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model,imaxRetries); 
% **
% **
% **	Input: 			mSigmaU: 	reduced form resdidual covariance matrix
% **				
% **					mS:			restriction matrix
% **	
% **					vs: 		vector of normalizing constants
% **
% **					vg: 		vector of start values
% **
% **					T:			number of obs. used in estimation
% **
% **					eps1_tol: 	tolerance for relative change in parameters
% **
% **					eps2_tol:	tolerance for relative change in log likelihood
% **
% **					DegreeofOverIdent: degree of overidentification
% **
% **					maxIterations:	max. number of iteration in scoring
% **
% **					model:      	scalar, model indicator
% **									1 = AB-model
% **									2 = A-model
% **									3 = b-model
% **
% **					imaxRetries:	max. # of redrawing new starting values for correlation matrix decomposition
% **
% **
% **		Output:		vgn:		optimal start values 
% **
% **					i_irt		indicator: 1 # of max. retries exceeded
% **
% **					i:			# of iterations needed 	
% **				
% **					noConvergence: 0 if convergence, 1 else
% **
% */
K = size(mSigmaU,1);
mD = zeros(K,K);
mD(1:K+1:end) = sqrt(diag(mSigmaU));
mCor = inv(mD)*mSigmaU*inv(mD);   	 
j  = 1;
noConvergence =1;
while noConvergence == 1 || j<=imaxRetries
    [mA,mB,~,~,~,~,~,i,~,~,~, noConvergence] = MLEstSVAR_svar_var(mCor,mS,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,model,0);  	   					
    if noConvergence == 1;
        vg = rndn(cols(mS),1);
    end
	j = j +1;
end

vgn = GetNewStartValues_svar_var(mA,mB,mS,mD,model);

if noConvergence == 1 && j>imaxRetries;
	i_irt = 1;
else
	i_irt = 0;
end	