function z=commutation_svar_var(n);		
%/* computes commutation matrix */
z = zeros(n^2,n);
c=1;			%/* first n columns */
r=1;
while c<=n
	z(r,c)=1;
    c=c+1;
    r=r+n; 
end       	   	
zold= z;
cc = 1;		  %     /* rest of colums  */
while cc<n
    z = [z [zeros(cc,n);zold(1:n^2-cc,:)]];
    cc = cc +1;
end
