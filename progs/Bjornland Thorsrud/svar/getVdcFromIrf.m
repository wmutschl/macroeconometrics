function vdc=getVdcFromIrf(irf)
% PURPOSE: Compute variance decompositions from impulse responses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input: 
%
% irf = array (n x n x h x d) of impulse resonses from the function
% varImpulseResponsFast. n is the number of variables in the VAR, h is the
% number of horizons, and d is the number of simulations (can be equal to
% 1)
%
% Output: 
%
% vdc = array, same size as irf. See function varImpulseResponsFast for
% further comments. 
%
% Usage: 
%
% vdc=getVdcFromIrf(irf)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nvary,nvars,hor,draws]=size(irf);
allirfs=zeros(nvary,1,draws);
sumirfs=zeros(nvary,nvars,draws);
vdc=nan(nvary,nvars,hor,draws);
for i=1:hor
    allirfs=allirfs+permute(sum(irf(:,:,i,:).^2,2),[1 2 4 3]);       
    sumirfs=sumirfs+permute(irf(:,:,i,:).^2,[1 2 4 3]);
    vdc(:,:,i,:)=sumirfs./repmat(allirfs,[1 nvars 1]);    
end;                                                


