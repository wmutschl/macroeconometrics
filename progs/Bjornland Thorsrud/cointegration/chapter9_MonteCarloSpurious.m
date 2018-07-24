function chapter9_MonteCarloSpurious()
%% PURPOSE: Monte Carlo illustration of spurious regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The code is not written for computational efficiency or
% elegance. The book: 
%
% "Applied Time Series for Macroeconomics"
% Gyldendal Akademisk 2014
% by Hilde C. Bjørnland and Leif A. Thorsrud 
%
% provides details. Please refer to the book if the code(s) are used for 
% research of commercial purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set the randstream (so that we get the same random numbers every time)
s=RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

%% Settings
N=5000;
T=[50 100 500];
xDomainLength=200;

%% Do the correct case 
xDomainb=linspace(-0.75,0.75,xDomainLength);
xDomaint=linspace(-10,10,xDomainLength);

% Empty output
b1iid=nan(xDomainLength,numel(T));
tvaliid=nan(xDomainLength,numel(T)+1);
dwd=nan(1,numel(T));

% Generate two iid N(0,1) variables
cnt=1;
for t=T    
    b1iidt=nan(N,1);
    tvaliidt=nan(N,1);
    dwdt=nan(N,1);
    for n=1:N                        
        
        x1=randn(t,1);
        x2=randn(t,1);
        
        const=ones(t,1);
        y=x1;
        x=[const x2];
        % regress: x1 = b0 +  b1*x2 + u 
        b0=x\y;
        % get tval
        resid=y-x*b0;        
        tmp=sqrt(diag((resid'*resid)/(t-2)));
        sige=tmp(:);                
        sigma=sige(:).^2;
        invx=diag(inv(x'*x));
        bstd=sqrt(sigma*invx);                                    
        tval=b0./bstd;            
        % allocate to output
        b1iidt(n,1)=b0(2);
        tvaliidt(n,1)=tval(2);                
       % [~,dwdt(n,1)]=dwtest(resid,x);
    end
    b1iid(:,cnt)=ksdensity(b1iidt,xDomainb);
    tvaliid(:,cnt)=ksdensity(tvaliidt,xDomaint);    
    dwd(1,cnt)=mean(dwdt);
    
    cnt=cnt+1;
end
tvaliid(:,cnt)=ksdensity(randn(N,1),xDomaint);

%% Output figures

figure('name','Beta convergence')
plot(xDomainb,b1iid(:,1),'k-','linewidth',2);
hold on
plot(xDomainb,b1iid(:,2),'k--','linewidth',2);
plot(xDomainb,b1iid(:,3),'color',[0.5 0.5 0.5],'linewidth',2);
lh=legend('T=50','T=100','T=500');
set(lh,'box','off');

figure('name','t-values')
plot(xDomaint,tvaliid(:,1),'k-','linewidth',2);
hold on
plot(xDomaint,tvaliid(:,2),'k--','linewidth',2);
plot(xDomaint,tvaliid(:,3),'k-*','linewidth',2);
plot(xDomaint,tvaliid(:,4),'color',[0.5 0.5 0.5],'linewidth',2,'marker','*');
lh=legend('T=50','T=100','T=500','N(0,1)');
set(lh,'box','off');
set(gca,'xlim',[-40 40]);

%% Do the spurious case 
xDomainb=linspace(-4,4,xDomainLength);
xDomaint=linspace(-75,75,xDomainLength);

% empty output
b1iidsc=nan(xDomainLength,numel(T));
tvaliidsc=nan(xDomainLength,numel(T));
dwdsc=nan(1,numel(T));

% Generate to iid N(0,1) variables
cnt=1;
for t=T    
    b1iidt=nan(N,1);
    tvaliidt=nan(N,1);
    dwdt=nan(N,1);
    for n=1:N                        
        x1=randn;
        x2=randn;
        for tt=2:t
            x1=cat(1,x1,0 + x1(tt-1)+randn);
            x2=cat(1,x2,0 + x2(tt-1)+randn);
        end
        
        const=ones(t,1);
        y=x1;
        x=[const x2];
        % regress: x1 = b0 +  b1*x2 + u 
        b0=x\y;
        % get tval
        resid=y-x*b0;        
        tmp=sqrt(diag((resid'*resid)/(t-2)));
        sige=tmp(:);                
        sigma=sige(:).^2;
        invx=diag(inv(x'*x));
        bstd=sqrt(sigma*invx);                                    
        tval=b0./bstd;            
        % allocate to output
        b1iidt(n,1)=b0(2);
        tvaliidt(n,1)=tval(2);                
        %[~,dwdt(n,1)]=dwtest(resid,x);
    end
    b1iidsc(:,cnt)=ksdensity(b1iidt,xDomainb);
    tvaliidsc(:,cnt)=ksdensity(tvaliidt,xDomaint);    
    dwdsc(1,cnt)=mean(dwdt);
    
    cnt=cnt+1;
end

%% Output figures

figure('name','Beta convergence - spurious')
plot(xDomainb,b1iidsc(:,1),'k-','linewidth',2);
hold on
plot(xDomainb,b1iidsc(:,2),'k--','linewidth',2);
plot(xDomainb,b1iidsc(:,3),'color',[0.5 0.5 0.5],'linewidth',2);
lh=legend('T=50','T=100','T=500');
legend(lh,'boxoff')

figure('name','t-values - spurious')
plot(xDomaint,tvaliidsc(:,1),'k-','linewidth',2);
hold on
plot(xDomaint,tvaliidsc(:,2),'k--','linewidth',2);
plot(xDomaint,tvaliidsc(:,3),'color',[0.5 0.5 0.5],'linewidth',2);
lh=legend('T=50','T=100','T=500');
legend(lh,'boxoff')
set(gca,'xlim',[-40 40]);


