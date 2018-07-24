function ystar = BootstrapGDP(VAR)
% =======================================================================
% Computes a standard residual-based bootstrap GDP for VAR models with
% iid sampling with replacement for errors and random sampling of block
% initial values
% =======================================================================
% ystar = BootstrapGDP(VAR)
% -----------------------------------------------------------------------
% INPUTS
%	- VAR     : Structure of reduced-form estimation function VARReducedForm.m. [structure]
% -----------------------------------------------------------------------
% OUTPUTS
%   - ystar   : Data matrix. [number of periods x number of variables]
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
% =======================================================================

ENDO = VAR.ENDO;
[T,K]=size(ENDO);
p = VAR.nlag;
U=VAR.residuals;
Acomp = VAR.Acomp;
vhat = zeros(size(Acomp,1),1);
dhat = zeros(size(Acomp,1),1);

if VAR.const == 1
    vhat(1:K,1) = VAR.A(:,1);
elseif VAR.const == 2
    vhat(1:K,1) = VAR.A(:,1);
    dhat(1:K,1) = VAR.A(:,2);
end

y=ENDO';
Y=y(:,p:T);	
for i=1:p-1
 	Y=[Y; y(:,p-i:T-i)];		
end

Ur=zeros(K*p,T-p);   
Yr=zeros(K*p,T-p+1); 

% iid resampling for blockwise initial values
pos=fix(rand(1,1)*(T-p+1))+1;
Yr(:,1)=Y(:,pos);

% iid resampling for error terms
index=fix(rand(1,T-p)*(T-p))+1;
Ur(:,2:T-p+1)=U(:,index);	

lintrend =  1:T;
for i=2:T-p+1
    Yr(:,i)= vhat + dhat.*lintrend(i) + Acomp*Yr(:,i-1)+Ur(:,i); 
end

ystar=[Yr(1:K,:)];
for i=2:p
    ystar=[Yr((i-1)*K+1:i*K,1) ystar];
end
ystar=ystar';
end