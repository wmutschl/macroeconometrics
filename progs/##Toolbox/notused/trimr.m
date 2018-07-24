%=========================================================================
% Returns a matrix (or vector) stripped of the specified rows
%
%   Inputs: 
%             x  = input matrix (or vector) (n x k)
%             rb = first n1 rows to strip
%             re = last  n2 rows to strip
%
%  Mimics the Gauss routine. 
%=========================================================================


function z = trimr(x,rb,re)


    n = length(x);
    if (rb+re) >= n; 
        error('Attempting to trim too much');
    end
    z = x(rb+1:n-re,:);
    
end
  
