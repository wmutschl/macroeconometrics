% FIGURE12_5.M
%
% Kilian and Lutkepohl (2017), Structural VAR Analysis, Cambridge University Press.
% This file generates Figure 12.5
clear; 

trivar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VAR Impulse response analysis 
[IRF]=irfvar(A,SIGMA(1:q,1:q),p,h);
    IRF(1,:)=cumsum(IRF(1,:));
    IRF(4,:)=cumsum(IRF(4,:));
    IRF(7,:)=cumsum(IRF(7,:));

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
    [Ar,SIGMAr]=olsvarc(yr,p);

    % Compute IRFs
    IRFr=irfvar(Ar,SIGMAr(1:q,1:q),p,h);
    IRFr(1,:)=cumsum(IRFr(1,:));
    IRFr(4,:)=cumsum(IRFr(4,:));
    IRFr(7,:)=cumsum(IRFr(7,:));
    IRFmat(r,:)=vec(IRFr)';
end;
IRFrstd=reshape((std(IRFmat)'),q^2,h+1);
CI1LO=IRF-1.96*IRFrstd; CI1UP=IRF+1.96*IRFrstd;

horizon=0:h;

subplot(3,3,1)
plot(horizon,-IRF(1,:),'r-',horizon,-CI1LO(1,:),'b:',horizon,-CI1UP(1,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Oil production','fontsize',16)
axis([0 h -2 0.5])

subplot(3,3,2)
plot(horizon,-IRF(2,:),'r-',horizon,-CI1LO(2,:),'b:',horizon,-CI1UP(2,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Real activity','fontsize',16)
axis([0 h -5 5])

subplot(3,3,3)
plot(horizon,-IRF(3,:),'r-',horizon,-CI1LO(3,:),'b:',horizon,-CI1UP(3,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil supply shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
axis([0 h -0.075 0.075])

subplot(3,3,4)
plot(horizon,IRF(4,:),'r-',horizon,CI1LO(4,:),'b:',horizon,CI1UP(4,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Oil production','fontsize',16)
axis([0 h -1.5 1])

subplot(3,3,5)
plot(horizon,IRF(5,:),'r-',horizon,CI1LO(5,:),'b:',horizon,CI1UP(5,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Real activity','fontsize',16)
axis([0 h -2 8])

subplot(3,3,6)
plot(horizon,IRF(6,:),'r-',horizon,CI1LO(6,:),'b:',horizon,CI1UP(6,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Aggregate demand shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
axis([0 h -0.05 0.1])

subplot(3,3,7)
plot(horizon,IRF(7,:),'r-',horizon,CI1LO(7,:),'b:',horizon,CI1UP(7,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Oil production','fontsize',16)
xlabel('Months','fontsize',16)
axis([0 h -1.5 1])

subplot(3,3,8)
plot(horizon,IRF(8,:),'r-',horizon,CI1LO(8,:),'b:',horizon,CI1UP(8,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Real activity','fontsize',16)
xlabel('Months','fontsize',16)
axis([0 h -5 5])

subplot(3,3,9)
plot(horizon,IRF(9,:),'r-',horizon,CI1LO(9,:),'b:',horizon,CI1UP(9,:),'b:',horizon,zeros(size(horizon)),'k-','linewidth',3);
title('Oil-specific demand shock','fontsize',16)
ylabel('Real price of oil','fontsize',16)
xlabel('Months','fontsize',16)
axis([0 h 0 0.15])
