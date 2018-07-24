%/* procedure that takes y, x and d matrix as input und returns 
%written by Markus Kraetzig 2001
%
%infocrits and a vector of optimal laglength 
%maxlag -- to search from  startlag  up to maxlag
%exlag -- lags of exogenous variables
%yr -- TxK matrix of endogenous variables
%x -- TxS matrix of exogenous variables
%d -- TxM matrix of deterministic variables

function [infocrits,crits] = var_infocrit(yr, x, const, startlag, maxlag,exlag)
aic_out = []; fpe_out=[]; hq_out=[]; sc_out=[];
if const == 0;
    d = [];
elseif const == 1
    d = ones(size(yr,1),1);
elseif const == 2
    d = [ones(size(yr,1),1) transpose(1:size(yr,1))];
elseif const == 3
    d = [ones(size(yr,1),1) transpose(1:size(yr,1)) transpose(1:size(yr,1)).^2];
end

%if x == 0;
%  x = [];
%end;

if exlag <= maxlag
    max = maxlag;
else
    if size(x,1) > 0
        max = exlag;
    else
        max = maxlag;
    end
end

y = yr(max+1:size(yr,1),:);
t = size(y,1);
j = startlag;
dime=size(yr,2);
while j <= maxlag
    %z = ones(1,t);
    if size(d,1) > 0
        z = d(max+1:size(d,1),:)';
    end
    if size(x,1) > 0
        i = 0;
        while i <= exlag
            z = [z x(max+1-i:size(x,1)-i,:)'];
            i = i+1;
        end
    end
    if size(y,1) > 0
        i = 0;
        while i < j
          z = [z;yr(max-i:size(yr,1)-1-i,:)'];
          i = i+1;
        end
    end
%     if j>0 || size(x,1)>0 || size(d,1)>0
%         z=transpose(trimr(z',1,0));
%     end
    b = inv(z*z')*z*y;
    resid = y'-b'*z;
    sigma=(resid*resid')/t;
    aic = log(det(sigma))+2*j*dime^2/t;
    aic_out= [aic_out aic];
    fpe = ((t+size(z,1))/(t-size(z,1)))^dime*det(sigma);
    fpe_out = [fpe_out fpe];
    sc = log(det(sigma))+log(t)/t*j*dime^2;
    sc_out=[sc_out sc];
    hq = log(det(sigma))+2*log(log(t))/t*j*dime^2;
    hq_out=[hq_out hq];
    j=j+1;
end

% aic_out=aic_out(2:rows(aic_out),.];
% fpe_out=fpe_out[2:rows(fpe_out),.];
% sc_out=sc_out[2:rows(sc_out),.];
% hq_out=hq_out[2:rows(hq_out),.];

alags=find(min(aic_out)) ,aic_out)-1;
flags=indnv(minc(fpe_out),fpe_out)-1;
hlags=indnv(minc(hq_out),hq_out)-1;
slags=indnv(minc(sc_out),sc_out)-1;
opt=alags|flags|hlags|slags;
opt = opt + startlag;
retp(aic_out~fpe_out~hq_out~sc_out,opt);
endp;


/* just give the values of the information criteria
input:
resids = T x K of resids from VAR
p = used lags ins VAR
n_of_para = the number of parameters estimated in the model
in the FPE, the average estimated coefficients in each equation are used
output: 
inf = AIC|FPE|SC|HQ
*/
proc(1)=display_infocrit(resids, n_of_para);
local inf,t,dime,sigma;
inf=zeros(4,1);
dime=cols(resids);
t=rows(resids);
sigma=(resids'*resids)./t;

inf[1,1] = ln(det(sigma))+2*n_of_para/t;

inf[2,1] = ((t+n_of_para/dime)/(t-n_of_para/dime))^dime*det(sigma);

inf[3,1] = ln(det(sigma))+ln(t)/t*n_of_para;

inf[4,1] = ln(det(sigma))+2*ln(ln(t))/t*n_of_para;

retp(inf);
endp;