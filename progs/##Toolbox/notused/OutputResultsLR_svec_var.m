function OutputResultsLR_svec_var(mB,mSigmaUt,i,LLnew,LR_stat,LR_prob,DegreeofOverIdent,T,noConvergence,maxIterations,vg,model,i_ident,iCorr,i_iRT,i_itercorr,i_coCorr,vg_corr,mC1ident,nResB, nResC1,BS_ind,mSE_mB, mt_mB,mSE_mC1ident,mt_mC1ident, fname)
% /*
% **	Purpose:	Write results of SVAR with long run restrictions to screen and ASCII-file
% **
% **		BS_ind: if 1, bootstrap std. errors and t-values are printed
% **				if 0, only point estimates are printed
% */
% proc(0) = OutputResultsLR_svec_var(mB,mSigmaUt,i,LLnew,LR_stat,LR_prob,DegreeofOverIdent,T,noConvergence,maxIterations,vg,model,i_ident,iCorr,i_iRT,i_itercorr,i_coCorr,vg_corr,mC1ident,nResB, nResC1,BS_ind,mSE_mB, mt_mB,mSE_mC1ident,mt_mC1ident, fname);
% local K;
% screen off;
% output file = ^fname reset;
% "";
% if model == 3;
% "This is a B-model with long run restrictions";
% endif;
% 
% format /rd 0,0;
% 
% if i_ident >= 0;
%  "";
%  "Long run restrictions provide(s) " nResC1 "independent restriction(s).";
%  "Contemporaneous restrictions provide(s) " nResB "additional restriction(s).";
%  "";
% endif;
% 
% K = rows(mC1ident);
% 
% format /rd 8,4;
% 
%   if i_ident == -1;
%     "";
%     "Invalid restrictions, implied B matrix singular!";
%   elseif i_ident == 0;  
% 
% 	if (nResC1+nResB) lt (K*(K-1)/2);
% 		"You have specified too few independent restrictions...";
% 		
% 	else;
% 	     "Model is not identified at starting values";
%   	endif;
%   	"";
%   	"Model is not identified, use different restrictions and try again!";
% else;
% 	if iCorr == 1 and i_iRT == 0 and i_coCorr == 0;
% 			"";
% 			"Step 1:";
% 			"Obtaining starting values from decomposition of correlation matrix...";
% 			"Iterations needed for correlation matrix decomposition: " i_itercorr;
% 			"Vector of rescaled starting values: " vg_corr;
% 			"";
% 			"Step 2:";
% 					if noConvergence == 0;
% 					format /rd 8,4;
% 					"Structural VAR Estimation Results";
% 					"ML Estimation, Scoring Algorithm (see Amisano & Giannini (1992))";
% 					
% 						format /rd 0,0; 
% 						"Convergence after " i "iterations";
% 						format /rd 8,4;
% 						"Log Likelihood: " LLnew;	
% 						if DegreeofOverIdent == 0;
% 						"Structural VAR is just identified";
% 						else;
% 						"Structural VAR is over-identified with " DegreeofOverIdent " degrees of freedom";
% 						"LR Test: Chi^2(" DegreeofOverIdent "): " LR_Stat ", Prob: " LR_prob ;
% 						endif;
% 						"";		
% 						"Estimated B matrix:"  mB;
% 						"";
% 						if BS_ind == 1;
% 							"Bootstrap standard errors:"
% 							mSE_mB;
% 							"";
% 							"Bootstrap t-values:"
% 							mt_mB;
% 							"";
% 						endif;	
% 						"Estimated long run impact matrix" mC1ident;
% 						"";
% 						if BS_ind == 1;
% 							"Bootstrap standard errors:"
% 							mSE_mC1ident;
% 							"";
% 							"Bootstrap t-values:"
% 							mt_mC1ident;
% 							"";
% 						endif;	
% 
% 						"SigmaU~*100" mSigmaUt*100;
% 						"end of ML estimation";
% 						"";
% 					else;
% 						"";
% 						"Warning!";
% 						format /rd 0,0;
% 						"No Convergence after "  maxIterations " iterations";
% 						format /rd 8,4;
% 						"Try again using different starting values...";
% 						"";
% 					endif;	
% 		elseif iCorr == 1 and (i_iRT == 1 or i_coCorr == 1);
% 		
% 			"Decomposition of correlation matrix failed";
% 			"Try again...";
% 		elseif iCorr == 0 and noConvergence == 0;
% 			format /rd 8,4;
% 			"Structural VAR Estimation Results";
% 			"ML Estimation, Scoring Algorithm (see Amisano & Giannini (1992))";
% 
% 				format /rd 0,0; 
% 				"Convergence after " i "iterations";
% 				format /rd 8,4;
% 				"Log Likelihood: " LLnew;	
% 				if DegreeofOverIdent == 0;
% 				"Structural VAR is just identified";
% 				else;
% 				"Structural VAR is over-identified with " DegreeofOverIdent " degrees of freedom";
% 				"LR Test: Chi^2(" DegreeofOverIdent "): " LR_Stat ", Prob: " LR_prob ;
% 				endif;
% 				"";		
% 				"Estimated B matrix"  mB;
% 				"";
% 				if BS_ind == 1;
% 						"Bootstrap standard errors:"
% 						mSE_mB;
% 						"";
% 						"Bootstrap t-values:"
% 						mt_mB;
% 						"";
% 				 endif;	
% 					"Estimated long run impact matrix" mC1ident;
% 					"";
% 					if BS_ind == 1;
% 						"Bootstrap standard errors:"
% 						mSE_mC1ident;
% 						"";
% 						"Bootstrap t-values:"
% 						mt_mC1ident;
% 						"";
% 					endif;	
% 				"SigmaU~*100" mSigmaUt*100;
% 				"end of ML estimation";
% 				"";
% 			elseif iCorr == 0 and noConvergence == 1;
% 				"";
% 				"Warning!";
% 				format /rd 0,0;
% 				"No Convergence after "  maxIterations " iterations";
% 				format /rd 8,4;
% 				"Try again using different starting values...";
% 				"";
% 			endif;
% endif;	/* endif i_ident*/
% output off;
% screen on;
% endp;

