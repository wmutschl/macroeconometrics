% This function calculates the vech of a sqauared matrix
function B = vech(A)

n      = size(A,1);
B      = zeros(n*(n+1)/2,1);
index  = 0;
for i=1:n
    B(1+index:index+n+1-i,1) = A(i:n,i);    
    index = index + n+1-i;    
end

