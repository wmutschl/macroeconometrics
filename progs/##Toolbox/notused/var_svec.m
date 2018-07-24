% * mSigmaU    - K x K matrix, reduced form covariance matrix
% * T          - number of observations used in estimation
% /**
%  * Restriction on C and PSI can be specified here
%  *
%  *   use:
%  *    '*' = unrestricted estimate 
%  *    '0'   = restricted to zero
%  *    'other numeric value' = restricted to specified value 

%mSigmaU = round(vecm.r1.EstCov,4);
%alpha = round(vecm.r1.paramVals.A*vecm.r1.paramVals.B(1),4); 
%beta = round(vecm.r1.paramVals.B./vecm.r1.paramVals.B(1),4);
%mGamma = round([vecm.r1.paramVals.B1 vecm.r1.paramVals.B2],4);
mSigmaU = VECM.EstCov;
alpha = VECM.paramVals.A;%*vecm.paramVals.B(1); 
beta = VECM.paramVals.B;%./vecm.paramVals.B(1);
mGamma = [VECM.paramVals.B1 VECM.paramVals.B2];

T = size(VECM.res,1);
mB_Res = repmat(-1e12,4,4); mB_Res(4,2) = 0;
mC1_Res = repmat(-1e12,4,4); mC1_Res(1,2:3) = 0; mC1_Res(:,4) = 0;

eps1_tol = 1E-6;	   %tolerance of relative change in parameters 
eps2_tol = 1E-10; 	   %tolerance for relative change in logLik
maxIterations = 500;
iStartValueMethod = 1; %1: draw randomly, 2: fixed values set by "fixStart", 3:  vector specicied by user "vStartUser"
fixStart = .1;         %fixed starting value, only needed when  iStartValueMethod  = 2 
vStartUser = [.1,.2,.4,.6,.7]'; %vector of start values, only needed when iStartValueMethod = 3  
%                                 dimension changes, rows(vStartUser) = cols(mS)	
iCorr = 0;             %decompose correlation matrix first to obtain SV	
imaxRetries	= 10;      %max. retries of randomly drawing starting values 

SVAR_UNRESTRICTED = -1e12;
K  = size(mSigmaU,1);
%bet=beta';

tmp=vec(mB_Res);
for i = 1:size(tmp,1)
    if tmp(i) < SVAR_UNRESTRICTED;
        tmp(i) = SVAR_UNRESTRICTED;
    end
end
mB_Res = reshape(tmp,K,K);

tmp=vec(mC1_Res);
for i = 1:size(tmp,1)
    if tmp(i) < SVAR_UNRESTRICTED;
        tmp(i)=SVAR_UNRESTRICTED;
    end
end
mC1_Res = reshape(tmp,K,K);

%/* initialize some variables */

mC1		  = ComputeC1_svec_var(alpha, beta, mGamma);
mS=[];
vs=[];
vg_b=[];
vg_corr 	  = 0;	
i_iRt		  = 0;
i_itercorr	  = 0;
i_coCorr	  = 0;
mB		  = 0;
mC		  = 0;
mSigmaAB	  = 0;
mSigmaUt	  = 0;
i		  = 0;
LLnew		  = 0;
LR_stat		  = 0;
LR_prob		  = 0;
DegreeofOverIdent = 0;
noConvergence	  = 0;
BS_int  = 0;

										
[mS,vs,model,nfreea,nfreeb,nfreeall, nResB, nResC1] = GetResMatricesLR_svec_var(mB_Res, mC1_Res, mC1);
if model < 0;
    return
end
vg  = StartingValues_svar_var(mS,iStartValueMethod,vStartUser,fixStart); 
[i_ident, DegreeofOverIdent] = CheckIdentificationLR_svec_var(vg,mS,vs,nfreeall, nResB, nResC1);

if i_ident == 1;     %/* only proceed when model is identified, otherwise inform user*/    		
    if iCorr == 1;	 %/* decompose correlation matrix first to obtain SV 	*/
        [vg_corr,i_iRT,i_itercorr, i_coCorr]  = DecompCorrelation_svar_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model, imaxRetries);	
        vg = vg_corr; 
        if i_iRT == 0 && i_coCorr == 0;	%/* only estimate when 1st step was successful*/  
            [mB,vg_b,mSigmaUt,i,LLnew,LR_stat,LR_prob, noConvergence] = MLEstSVECM_svec_var(mSigmaU, mS,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,model);
        end
    else
        [mB,vg_b,mSigmaUt,i,LLnew,LR_stat,LR_prob, noConvergence] = MLEstSVECM_svec_var(mSigmaU, mS,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,model);
    end
end

mCB_point = mC1*mB;
mB_point  = mB;
vg_BS 	  = vg_b;

%/* Std. deviations and t-ratios */

m_se_B=[];
m_tv_B=[];
m_se_mCB=[];
m_tv_mCB=[];

%/* Bootstrap std deviations */
if BS_int  
  [m_se_B,m_se_mCB,m_tv_B,m_tv_mCB] = BootstrapStdErr_svec_var(var,beta,beta_d,bootRep,seed,mB_Res,mC1_Res,vg_bs,T,K,mB_point,mCB_point,eps1_tol, eps2_tol,DegreeofOverIdent, maxIterations);	
end

%OutputResultsLR_svec_var(mB,mSigmaUt,i,LLnew,LR_stat,LR_prob,DegreeofOverIdent,T,noConvergence, maxIterations,vg,model,i_ident,iCorr, i_iRT,i_itercorr, i_coCorr,vg_corr,mCB_point, nResB, nResC1,BS_int,m_se_B,m_tv_B,m_se_mCB,m_tv_mCB,fname);