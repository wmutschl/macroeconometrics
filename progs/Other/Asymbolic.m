K = 4;
p = 2;
const = 1;

c = sym('c',[K,1]);
d = sym('d',[K,1]);
if const == 0
    A = [];
elseif const == 1
    A = c;
elseif const == 2
    A = [c d];
end
for i = 1:p
    A = [A sym(sprintf('A%d_',i),[K,K])];
end
a = A(:)
V_diag = diag(diag(a*transpose(a)))