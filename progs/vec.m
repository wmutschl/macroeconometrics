% Willi Mutschler
function x=vec(y)
x=[];
for i=1:cols(y)
    x=[x;y(:,i)];
end