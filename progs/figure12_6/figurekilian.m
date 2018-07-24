% FIGUREKILIAN.M

clear; 

global h t

% Load monthly data for 1973.2-2007.12 ordered as:
% 1. Growth rate of world oil production 
% 2. Global real activity (index based on dry cargo shipping rates)
% 3. Real price of oil 
% The data sources are described in the text.

load data.txt; y=data; [t,q]=size(y); 
time=(1973+2/12:1/12:2007+12/12)';  % Time line
h=15;                               % Impulse response horizon
p=24;                               % VAR lag order
y=detrend(y,0);
[A,SIGMA,Uhat,V,X]=olsvarc(y,p);	% VAR with intercept	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR Impulse response analysis 
[IRF]=irfvar(A,SIGMA(1:q,1:q),p,h);
    IRF(1,:)=cumsum(IRF(1,:));
    IRF(4,:)=cumsum(IRF(4,:));
    IRF(7,:)=cumsum(IRF(7,:));
    
if ~ any(abs(eig(A))>=1)
	[A]=asybc(A,SIGMA,t,p);
end;
SIGMA=SIGMA(1:q,1:q);

% VAR bootstrap
randn('seed',1234);
nrep=2000;
IRFmat=zeros(nrep,(q^2)*(h+1));

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
    r    
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
    yr=detrend(yr,0);

    [Ar,SIGMAr]=olsvarc(yr,p);
    if ~ any(abs(eig(Ar))>=1)
		[Ar]=asybc(Ar,SIGMAr,t,p);
	end;


    % Compute IRFs
    IRFr=irfvar(Ar,SIGMAr(1:q,1:q),p,h);
    IRFr(1,:)=cumsum(IRFr(1,:));
    IRFr(4,:)=cumsum(IRFr(4,:));
    IRFr(7,:)=cumsum(IRFr(7,:));
    IRFmat(r,:)=vec(IRFr)';
end;

% Calculate 68/95 percent interval endpoints
%CI=prctile(IRFmat,[16 84]);
CI=prctile(IRFmat,[2.5 97.5]);
CILO=reshape(CI(1,:),q^2,h+1);
CIUP=reshape(CI(2,:),q^2,h+1);


horizon=0:h;

subplot(3,3,1)
plot(horizon,-IRF(1,:),'r-',horizon,-CILO(1,:),'b--',horizon,-CIUP(1,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Oil production','fontsize',16)
axis([0 h -2 0.5])

subplot(3,3,2)
plot(horizon,-IRF(2,:),'r-',horizon,-CILO(2,:),'b--',horizon,-CIUP(2,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Real activity','fontsize',16)
axis([0 h -5 5])

subplot(3,3,3)
plot(horizon,-IRF(3,:),'r-',horizon,-CILO(3,:),'b--',horizon,-CIUP(3,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
axis([0 h -0.05 0.05])

subplot(3,3,4)
plot(horizon,IRF(4,:),'r-',horizon,CILO(4,:),'b--',horizon,CIUP(4,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Oil production','fontsize',16)
axis([0 h -1.5 1])

subplot(3,3,5)
plot(horizon,IRF(5,:),'r-',horizon,CILO(5,:),'b--',horizon,CIUP(5,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Real activity','fontsize',16)
axis([0 h -2 8])

subplot(3,3,6)
plot(horizon,IRF(6,:),'r-',horizon,CILO(6,:),'b--',horizon,CIUP(6,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
axis([0 h -0.02 0.08])

subplot(3,3,7)
plot(horizon,IRF(7,:),'r-',horizon,CILO(7,:),'b--',horizon,CIUP(7,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Oil production','fontsize',16)
xlabel('Months')
axis([0 h -1.5 1])

subplot(3,3,8)
plot(horizon,IRF(8,:),'r-',horizon,CILO(8,:),'b--',horizon,CIUP(8,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Real activity','fontsize',16)
xlabel('Months')
axis([0 h -5 5])

subplot(3,3,9)
plot(horizon,IRF(9,:),'r-',horizon,CILO(9,:),'b--',horizon,CIUP(9,:),'b--',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
xlabel('Months')
axis([0 h 0 0.1])

save kilian95
