function bicVal=bicTest(epe1,neqs,lags,nobst)
% PURPOSE: Comput BIC test evaluation criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bicVal=log(det(epe1))+(neqs*neqs*lags)*(log(nobst)/nobst);         
