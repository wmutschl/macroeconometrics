%=========================================================================
%
%   Program to mimic the GAUSS seqa command
%
%========================================================================
    
function  y = seqa(start,step,nvals)

    y = start:step:start+step*(nvals - 1);

end

