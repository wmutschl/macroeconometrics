
function y = recserar(x,y0,a)
%Computes a vector of autoregressive recursive series (Gauss translation)

[N,K]=size(x);
[P,k]=size(y0); %a is P*K
y=y0;
for t=[P+1:1:N];
   yhere=x(t,:);
   for i=1:P; yhere=yhere+a(i,:).*y(t-i,:); end;
   y=[y;yhere];
end;
