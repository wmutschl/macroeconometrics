% VERSION1.M
% Normalization: diag(B0)=ones(q,1) and SIGMA_w diagonal

function q=version1(guess)

global A SIGMA p

K=size(guess,1);
B0=guess(1:K,1:K); % Initial guess for B0
B0inv=inv(B0);

F=vec(B0inv*diag(guess(:,K+1),0)*B0inv'-SIGMA(1:K,1:K));

% Using recursive short-run restrictions (like in chol-decomposition)
q=[F;  B0inv(1,2);  B0inv(1,3);  B0inv(2,3); B0(1,1)-1; B0(2,2)-1;  B0(3,3)-1];


q'+1


