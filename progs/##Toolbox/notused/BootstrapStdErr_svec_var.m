function [m_se_B,m_se_mCB,m_tv_B,m_tv_mCB] = BootstrapStdErr_svec_var(var,beta,beta_d,bootRep,seed,mB_Res,mC1_Res,vg_bs,T,K,mB,mCB,eps1_tol, eps2_tol,DegreeofOverIdent, maxIterations);


% /*
% **	Computes bootstrap standard errors and t-values of
% **	contemp. and long run impact matrix. It also computes the point estimates again.
% **	
% **		      	
% **                      var:            databuffer with parameters of reduced form estimation
% **                      bootRep:        skalar, number of bootstrap replication 
% **                      seed:           skalar, seed value
% **
% **
% **	Output: 
% **
% **                      m_se_B:	        K x K, bootstrap std errors of contemp. impact matrix	
% **
% **			m_se_mCB 	K x K, bootstrap std errors of long run impact matrix	
% **
% **                      m_tv_B:	        K x K, bootstrap t-values of contemp. impact matrix	
% **
% **			m_tv_mCB 	K x K, bootstrap t-values of long run impact matrix	
% 
% **
% **			
% **	
% **
% */

m_se_B=[];
m_se_mCB=[];
m_tv_B=[];
m_tv_mCB=[];
%@ model estimation without estimating beta again @  
var = vml_vdel(var, "em"$|"r_est"$|"cir"$|"idx_equa");
var = vml_SetCointRelation(var, beta, beta_d);
varHat = var_EstimateModel(var);   

@ prepare all parameters to recompute the time series @
{A0,A,B,C_VAR,F_VAR, mx_c_var,mx_tf, y0,u,x,d,z} = 
vml_residualBootstrap_prepare(varHat);

mmB0  = zeros(K^2,bootRep);
mmC0  = zeros(K^2,bootRep);
del_boot_rep=zeros(bootRep,1);
  
	for i(1,bootRep,1);
	  @ compute bootstrap time series @
	  {y_star,seed} = 
	  vml_residualBootstrap(A0,A,B,C_VAR,F_VAR, mx_c_var,mx_tf, y0,u,x,d,z, seed);

	  @ update var @
	  varHat = vml_SetEndogenousVariables(var, y_star);
	  
	  @ estimate model (bootstrap data) @
	  varHat = var_EstimateModel(varHat);
	  	  	  
	  mSigmaU_BS    = vml_VeRead(varHat,"cvRes");
	  alpha_BS      = vml_VeRead(varHat,"alpha");
	  mGamma_BS     = vml_VeRead(varHat,"G");
	  
	  bet=beta';
	  mC1_BS = ComputeC1_svec_var(alpha_BS,bet,mGamma_BS);
	  
	  {mS_BS,vs_BS,model,nfreea,nfreeb,nfreeall, nResB, nResC1}  =  
	  GetResMatricesLR_svec_var(mB_Res, mC1_Res, mC1_BS); 	
	  
	  if cols(mS_BS) == rows(vg_BS) and model > -10;
	    mB_BS = GetStructDecompinBS_svec_var(mSigmaU_BS, mS_BS, vs_BS, vg_BS, T, eps1_tol, eps2_tol, 
	    DegreeofOverIdent, maxIterations,model,i);
	  
	    mCident_BS = mC1_BS*mB_BS;
	    
	    tmp = vec(mB_BS);
   	    if rows(mmB0[.,i]) == rows(tmp);	      
	      mmB0[.,i] = tmp;       /* store BS results in matrix*/
	    endif;  
	    tmp = vec(mCident_BS);
	    if rows(mmC0[.,i]) == rows(tmp);
	      mmC0[.,i] = tmp;  /* store BS results in matrix*/
	    endif;
	    
	    @ Progress notification @
	    var__ProgressNotification(i+0,bootRep);
	  else;
	    "failure in SVEC bootstrap replication" i;
	    del_boot_rep[i,1]=1;
	  endif;
	  
	endfor;
	
	/* remove failed bootstrap replications */
	if sumc(del_boot_rep) > 0;
	  mmB0=delif(mmB0', del_boot_rep)';
	  mmC0=delif(mmC0', del_boot_rep)';
	endif;
	
	/* compute bootstrap standard errors of B 	*/
	m_se_B     = reshape(sqrt(meanc(((mmB0-vec(mB))^2)')),K,K)';		
	
	/* compute bootstrap standard errors of C1 	*/
	m_se_mCB = reshape(sqrt(meanc(((mmC0-vec(mCB))^2)')),K,K)';	
	
	
	/* compute t-values and take care of zero division */
	m_tv_B=zeros(rows(mB),cols(mB));
	m_tv_mCB=zeros(rows(mB),cols(mB));
	for i(1,rows(mB),1);
	  for j(1,cols(mB),1);
	    if mB_res[i,j] ne 0;
	      m_tv_B[i,j]= mB[i,j]./m_se_B[i,j] ; 	  	 
	    endif;
	    if mC1_res[i,j] ne 0;	      
	      m_tv_mCB[i,j]= mCB[i,j]./m_se_mCB[i,j] ; 	  	 
	    endif;
	  endfor;
	endfor;