% VERSION2.M
% Normalization: SIGMAw=I

function q=version2(guess)

global A SIGMA p

K=size(guess,1);
B0inv=guess;

F=vec(B0inv*B0inv'-SIGMA(1:K,1:K));

% Using recursive short-run restrictions (like in chol-decomposition)
q=[F;  B0inv(1,2);  B0inv(1,3);  B0inv(2,3)];

q'+1


