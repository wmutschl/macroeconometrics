% BOOTSE.M

function IRFrse=bootse(y,p)

[t,q]=size(y); h=15;                        
[A,SIGMA,Uhat,V,X]=olsvarc(y,p); SIGMA=SIGMA(1:q,1:q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR bootstrap
randn('seed',1);
nrep=500;
IRFrmat=zeros(nrep,(q^2)*(h+1));

[t,q]=size(y);				
y=y';
Y=y(:,p:t);	
for i=1:p-1
 	Y=[Y; y(:,p-i:t-i)];		
end;

Ur=zeros(q*p,t-p);   
Yr=zeros(q*p,t-p+1); 
U=Uhat;    

for r=1:nrep
        
	pos=fix(rand(1,1)*(t-p+1))+1;
	Yr(:,1)=Y(:,pos);

    % iid resampling
	index=fix(rand(1,t-p)*(t-p))+1;
	Ur(:,2:t-p+1)=U(:,index);	

    for i=2:t-p+1
		Yr(:,i)= V + A*Yr(:,i-1)+Ur(:,i); 
	end;

	yr=[Yr(1:q,:)];
	for i=2:p
		yr=[Yr((i-1)*q+1:i*q,1) yr];
    end;
    yr=yr';
    [Ar,SIGMAr]=olsvarc(yr,p);
    
    % Compute IRFs
    IRFr=irfvar(Ar,SIGMAr(1:q,1:q),p,h);
    IRFr(1,:)=cumsum(IRFr(1,:));
    IRFr(4,:)=cumsum(IRFr(4,:));
    IRFr(7,:)=cumsum(IRFr(7,:));

    % Store results
    IRFrmat(r,:)=vec(IRFr)';
end;

% Compute standard deviation across r
IRFrse=std(IRFrmat);
