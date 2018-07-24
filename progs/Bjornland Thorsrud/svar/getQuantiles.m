function sortedValue=getQuantiles(value,alphaSign)
% PURPOSE: Get desired quantiles from simulated dirstribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
column=4;
draws=size(value,column);

selections=[ceil((alphaSign/2)*draws) round(0.5*draws) floor((1-alphaSign/2)*draws)];              
[sortedValue,order]=sort(value,column); 
sortedValue=sortedValue(:,:,:,selections);            