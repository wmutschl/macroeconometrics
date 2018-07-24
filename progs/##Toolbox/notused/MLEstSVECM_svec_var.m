function [mB,vg_B,mSigmaUt,i,llnew,LR_stat,LR_prob,noConvergence] = MLEstSVECM_svec_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model)
% /*
% **	Purpose: Estimates structural decomposition by ML based in VECM
% **	
% **  Usage:	{mB,mSigmaUt,i,LLnew,LR_stat,LR_prob,noConvergence}
% **			  =	MLEstSVAR(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model);
% **	
% **
% **	Input: 		
% **			mSigmaU		reduced form covariance matrix
% **
% **				mS			restriction matrix
% **
% **				vg			vector with starting values
% **
% **				T			observation used in estimation
% **
% **				eps1_tol	tolerance for relative parameter change
% **
% **				eps2_tol 	tolerance for relative change in log likelihood
% **
% **				DegreeofOverIdent		degreeofoveridentification
% **
% **				maxIterations	max. # of iterations
% **			
% **				model		model indicator: 1=AB, 2=K, 3=C
% **
% **	Output:						
% **				mB			estimated B matrix
% **
% **				vg_B		vector of free parameters, used as start values in bootstrap
% **
% **				mSigmaUt	ML estimate of covariance matrix
% **
% **				i 			# of iterations needed
% **
% **				LLnew		max LogLikelihood value
% **
% **				LR_stat		LR test statistic for overidentifying restrictions
% **				
% **				LR_prob		marginal significance level for LR_Stat
% **
% **			noConvergence	0 if alg. converged, else 1
% **
% */
maxls   = 1;
K 	 	= size(mSigmaU,1);
vecAB   = mS*vg+vs;
mA   	= reshape(vecAB(1:K^2,1),K,K);
mB   	= reshape(vecAB(K^2+1:2*K^2,1),K,K);
mCom 	= commutation_svar_var(K);
noConvergence = 0;

vgold   = vg;

eps1 	= 100;
eps2 	= 100;

mSigmaUt = inv(mA)*mB*mB'*inv(mA');
mK = inv(mB)*mA;
llold = T/2*log(det(mK)^2)-T/2*sum(diag((mK'*mK*mSigmaUt)));

i = 0;
while (isempty(find(eps1 > eps1_tol)) || (eps2 > eps2_tol)) && i<maxIterations		
    mK = inv(mB)*mA;
    mIAB = T*[ kron(inv(mK),inv(mB')) ; -1*kron(eye(K),inv(mB')) ]*(eye(K^2)+ mCom)*[ kron(inv(mK'),inv(mB)) -1*kron(eye(K),inv(mB))];
    mIga = mS'*mIAB*mS;
    v_scoreK 	 = T*(vec(inv(mK)')'-vec(mK)'*kron(mSigmaU,eye(K)));
    v_scoreAB 	 = v_scoreK * [kron(eye(K),inv(mB))  -1*kron(mA'*inv(mB'),inv(mB))];
    v_scoregamma = v_scoreAB*mS;  
   		   
   %/* RB, MK changed this block on May 22, 2003 */
   tmp ={};
%   trap 1;
   tmp = inv(mIga);
%   trap 0;
%   if scalerr(tmp);
%    retp(0,0,0,0,0,0,1);	/* return noConvergence = 1*/
%   endif;
   length = max(abs(tmp*v_scoregamma'));
   if length > maxls
    lambda = maxls/length;
   else
     lambda =1;
   end   

   vg = vgold +lambda*tmp*v_scoregamma';
	
    vecAB = mS*vg+vs;
	   	   
    mA = reshape(vecAB(1:K^2),K,K);
    mB = reshape(vecAB(K^2+1:2*K^2),K,K);   
    mK = inv(mB)*mA;
    mSigmaUt = inv(mA)*mB*mB'*inv(mA');
   
    llnew  = T/2*log(det(mK)^2)-T/2*sum(diag((mK'*mK*mSigmaUt)));

    eps2 = abs((llnew-llold)/llold);
	eps1 = abs((vg-vgold)./vgold);
    vgold = vg;	
    llold  = llnew;
    i = i +1;  
end
	

if isempty(find(eps1>eps1_tol))==0 || (eps2>eps2_tol);
  noConvergence = 1;
end
	
if sum(diag(inv(mA)*mB)<0) ~= 0;               %/* normalize sign of mB */
    x = diag(inv(mA)*mB)<0;
    mB(:,find(x==1)) = -1*mB(:,find(x==1));	
end

vg_B = inv(mS'*mS)*mS'*vec([mA mB]);

%/* compute LR test for over-identifying restr. */				    
if DegreeofOverIdent > 0;       
    LR_stat =  T*(log(det(mSigmaUt))-log(det(mSigmaU)));
    LR_prob =  chi2cdf(LR_stat,DegreeofOverIdent);
else
    LR_stat = 0;
    LR_prob = 0;
end
