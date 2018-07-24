proc(12)=MLEstSVAR_svar_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model,bs);
% /*
% **	Purpose: Estimates structural decomposition by ML
% **	
% **  Usage:	{mA,mB,minvAB,mSigmaAB,mSigmaUt,mA_std,mB_std,i,LLnew,LR_stat,LR_prob,noConvergence}
% **			  =	MLEstSVAR_svar_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model);
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
% **                              bs              1 when used by bootstrap (computes less parameters), 
% **                                              0 otherwise
% **
% **	Output:		mA			estimated A matrix
% **				
% **				mB			estimated B matrix
% **
% **				minvAB      inv(A)*B, contemp. impact matrix
% **
% **
% **				mSigmaAB	covariance matrix of structural parameter vector 
% **
% **				mSigmaUt	ML estimate of covariance matrix
% **
% **				mA_std		matrix with standard errors of mA
% **
% **				mB_std		matrix with standard errors of mB
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
K 	= size(mSigmaU,1);
vecAB   = mS*vg+vs;
mA   	= reshape(vecAB(1:K^2),K,K)';
mB   	= reshape(vecAB(K^2+1:2*K^2),K,K)';
mCom 	= commutation_svar_var(K);
noConvergence = 0;

vgold   = vg;

eps1 	= 100;
eps2 	= 100;

mSigmaUt = inv(mA)*mB*mB'*inv(mA');
   	  mK = inv(mB)*mA;
   llold = T/2*log(det(mK)^2)-T/2*sum(diag((mK'*mK*mSigmaUt)));

i = 0;
while (isempty(find(eps1>eps1_tol)) || (eps2>eps2_tol)) && i<maxIterations
		mK = inv(mB)*mA;            
        mIAB = T*(inv(mK).*.inv(mB')|-1*(eye(K).*.inv(mB')))*(eye(K^2)+ mCom)*(inv(mK').*.inv(mB)~-1*(eye(K).*.inv(mB)));        
        mIga = mS'*mIAB*mS;     
        v_scoreK 	 = T*(vec(inv(mK)')'-vec(mK)'*(mSigmaU.*.eye(K)));
        v_scoreAB 	 = v_scoreK * ((eye(K).*.inv(mB))~-1*mA'inv(mB').*.inv(mB));
        v_scoregamma = v_scoreAB*mS;  
   		   
   	   %/* RB, MK changed this block on May 22, 2003 */
   	   tmp ={};
   	   trap 1;
   	   tmp = inv(mIga);
   	   trap 0;
   	   if scalerr(tmp);
   	   	retp(0,0,0,0,0,0,0,0,0,0,0,1);	/* return noConvergence = 1*/
   	   endif;
       length = maxc(abs(tmp*v_scoregamma'));
       if length > maxls;
         lambda = maxls/length;
       else;
       	 lambda =1;
       endif;   
       
	  	 vg = vgold +lambda*tmp*v_scoregamma';
	
	   vecAB = mS*vg+vs;
	   
	   
 	    mA = reshape(vecAB[1:K^2],K,K)';
	    mB = reshape(vecAB[K^2+1:2*K^2],K,K)';	 
   
	   mK = inv(mB)*mA;
 mSigmaUt = inv(mA)*mB*mB'*inv(mA');
   
   LLnew  = T/2*ln(det(mK)^2)-T/2*sumc(diag((mK'mK*mSigmaUt)));	
   	 
     eps2 = abs((llnew-llold)/llold);
	 eps1 = abs((vg-vgold)/vgold);
    vgold = vg;	
   llold  = llnew;
   i = i +1;  
   
  endo;

if (eps1.>eps1_tol) or (eps2>eps2_tol);
  noConvergence = 1;
endif;  
	
	
	if model == 3;                                   	 /* for b-model          */
		if sumc(diag(inv(mA)*mB).<0) ne 0;               /* normalize sign of mB */
		  x = diag(inv(mA)*mB).<0; 					
		mB[.,indexcat(x,1)] = -1*mB[.,indexcat(x,1)];
		endif; 
    elseif model == 2;									 /* for A-model			 */
		if sumc(diag(inv(mA)*mB).<0) ne 0;               /* normalize sign of mA */
				  x = diag(inv(mA)*mB).<0; 		
				mA[.,indexcat(x,1)] = -1*mA[.,indexcat(x,1)];
		endif; 
		mA = NormSignK_svar_var(mA,mSigmaU);
	elseif model == 1;							/* for AB-model         */
	   if sumc(diag(mA).<0) ne 0;               /* normalize sign of mA */
		  x = diag(mA).<0; 				
		 mA[.,indexcat(x,1)] = -1*mA[.,indexcat(x,1)];
	 	endif; 
	
		if sumc(diag(mB).<0) ne 0;               /* normalize sign of mB */
			  x = diag(mB).<0; 				
			 mB[.,indexcat(x,1)] = -1*mB[.,indexcat(x,1)];
		endif;  
	
	      endif;
	      
	/* Quit here when in bootstrap */      
	if bs;
	  retp(mA,mB,inv(mA)*mB,0,0,0,0,i,0,0,0,0);
	endif;
	  
	/* compute LR test for over-identifying restr. */				    
      if DegreeofOverIdent > 0;       
        LR_stat =  T*(ln(det(mSigmaUt))-ln(det(mSigmaU)));
        LR_prob =  cdfchic(LR_stat,DegreeofOverIdent);
      else;
        LR_Stat = 0;
        LR_Prob = 0;
     endif;
	
	   /* compute standard errors */
           mK = inv(mB)*mA;            
   	 	 mIAB = T*(inv(mK).*.inv(mB')|-1*(eye(K).*.inv(mB')))*(eye(K^2)+ mCom)*(inv(mK').*.inv(mB)~-1*(eye(K).*.inv(mB)));        
 	     mIga = mS'*mIAB*mS; 
	mSigmaAB  = (mS*inv(mIga)*mS');
	std_AB    = sqrt(diag(mSigmaAB));
	   mA_std = reshape(std_AB[1:K^2],K,K)';
	   mB_std = reshape(std_AB[K^2+1:rows(std_AB)],K,K)';	 

retp(mA,mB,inv(mA)*mB,mSigmaAB,mSigmaUt,mA_std,mB_std,i,LLnew,LR_stat,LR_prob,noConvergence);
endp;
