function [ident,DegreeofOverIdent] = CheckIdentificationLR_svec_var(vg,mS,vs,nfreeall, nResB, nResC1)
% /* 
% **	Purpose:	Checks identification numerically using the starting values
% **
% **	Usage:		{ident, DegreeofOverIdent} = CheckIdentification(vg,mS,vs,nfreeall); 
% **
% **	Input:		vg: 		vector of starting values
% **				
% **				mS: 		Restriction matrix
% **				
% **				vs: 		vector of normalizing constants
% **	
% **				nfreeall:	scalar, number of free parameters
% **
% **	Output:		ident:		1: model identified, 0 model not identified
% **				
% **				DegreeofOverIdent:	gives the degree of overidentification
% */
ident = 0;
K = sqrt(size(vs,1)/2);
if (nResB+nResC1) < (K*(K-1)/2)   	    
    DegreeofOverIdent = 0;
else 	 
    vAB = mS*vg+vs;
	mA = reshape(vAB(1:K^2),K,K)';
	mB = reshape(vAB(K^2+1:2*K^2),K,K)';
	
	%/* build in trap for singularity of mB */
	   
% 	     tmp ={};
% 	     trap 1;
% 		 tmp = inv(mB);
%      	 trap 0;
% 	     if scalerr(tmp);
%       	       	/*errorlog("Invalid restrictions, implied B matrix singular!");*/
% 	     retp(-1,0);	/* return noConvergence = 1*/
%    	   	 endif;
	   
    mCom = commutation_svar_var(K); %commutation(K,K)
    %/* see page 173, Lütkepohl Kraetzig Applied TSA*/   
    mV   = [ (eye(K^2)+mCom)*kron((inv(mA)*mB)',inv(mB))  -1*(eye(K^2)+mCom)*kron(eye(K),inv(mB)) ]; 
    %mV   = (eye(K^2)+mCom)*((inv(mA)*mB)'.*.inv(mB)) ~ -1*(eye(K^2)+mCom)*(eye(K).*.inv(mB)); 
   	DegreeofOverIdent = K*(K+1)/2 - nfreeall;
	   
    if isempty(mV)==0 && isempty(mS)==0
        if sum(eig((mV*mS)'*(mV*mS)) < 1E-10) == 0
            ident = 1;
        end
    end	  	
end