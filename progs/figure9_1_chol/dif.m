% DIF.M
% Lutz Kilian
% University of Michigan
% April 1997
%
% This function returns the first differences of a t x q matrix of data

function [ydif]=dif(y)

t=size(y,1);
ydif=y(2:t,:)-y(1:t-1,:);

