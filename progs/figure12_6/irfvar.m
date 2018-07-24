% IRFVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

function [IRF]=irfvar(A,SIGMA,p,h)

J=[eye(3,3) zeros(3,3*(p-1))];
IRF=reshape(J*A^0*J'*chol(SIGMA)',3^2,1);

for i=1:h
	IRF=([IRF reshape(J*A^i*J'*chol(SIGMA)',3^2,1)]);
end;
