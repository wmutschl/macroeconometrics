function aicVal=aicTest(epe1,neqs,lags,nobst)
% PURPOSE: Comput AIC test evaluation criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aicVal=log(det(epe1))+(2*lags*neqs*neqs)/nobst;
