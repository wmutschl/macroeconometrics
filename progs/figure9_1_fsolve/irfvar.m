% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

function [IRF]=irfvar(A,B0inv,p,h)

q=size(B0inv,1);
J=[eye(q,q) zeros(q,q*(p-1))];
IRF=reshape(J*A^0*J'*B0inv,q^2,1);

for i=1:h
	IRF=([IRF reshape(J*A^i*J'*B0inv,q^2,1)]);
end;

