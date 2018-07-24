function Acomp = CompanionForm(p,A)
K = size(A,1);
Acomp = [A(:,2:end);
         eye(K*(p-1)) zeros(K*(p-1),p)];

