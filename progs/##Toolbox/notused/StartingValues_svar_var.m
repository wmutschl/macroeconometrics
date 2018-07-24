function vg = StartingValues_svar_var(mS,iStartValueMethod,vStartUser,fixStart)
% /*
% **	Purpose: Returns vector of start values with approriate dimension
% **
% **	Usage:		vg = StartingValues(mS,iStartValueMethod,vStartUser,fixStart); 
% **
% **	Input:		mS:	 restriction matrix, cols provide # of start values needed
% **
% **				iStartValueMethod:	1 - draw start values randomly
% **									2 - use fix start value for all parameters, e.g. .1
% **										contained in fixStart
% **									3 - use user specified vector of startvalues contained in 
% **										vStartUser
% **
% **				vStartUser:	vector of start values specified by user
% **
% **				fixStart:	scalar value of fixed starting value
% **
% **	Output:		vg - vector of starting values
% **
% */

if iStartValueMethod == 1
    vg = randn(size(mS,2),1)*.1;
elseif iStartValueMethod == 2
    vg = ones(size(mS,2),1)*fixStart;
elseif iStartValueMethod == 3
    vg = vStartUser;
end
