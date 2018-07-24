function chapter7_simulateVAR()
% PURPOSE: Simulate VAR data and compute some statistics
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
s=RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);

%% Settings VAR(1)
nlag=1;
nvary=2;
T=101;
sigma2=[1 0.2;0.2 1];
mu=[1;1];
phi=[0.5 0;1 0.2];
y0=[2;1];
% random erros (draw from multivariate normal using chol decomp)
e=chol(sigma2,'lower')*randn(nvary,T);

%% Simulate

y=[y0 nan(nvary,T-1)];
% generate VAR
for t=2:T    
    y(:,t)=mu+phi*y(:,t-1)+e(:,t);
end;
Y=y';

%% Compute some extra stuff

% Stability of var
if any(abs(real(eig(phi)))>=1)
    error('simulateVAR:stability','The VAR is not stationary given the parameters you have specified!')
end;

% Steady state of VAR
As=0;
for j=1:nlag
    As=As+phi(:,(j-1)*nvary+1:j*nvary);
end
SS=(eye(nvary)-As)\mu;                    

% Compute acf function
nacf=21;
[R,Gamma]=deal(zeros(nvary,nvary,nacf));   
Gamma(:,:,1)=reshape(inv(eye(nvary^2)-kron(phi,phi))*sigma2(:),[nvary nvary]);
D=diag(sqrt(diag(Gamma(:,:,1))));
R(:,:,1)=inv(D)*Gamma(:,:,1)*inv(D);
for i=2:nacf
    Gamma(:,:,i)=phi*Gamma(:,:,i-1);
    R(:,:,i)=inv(D)*Gamma(:,:,i)*inv(D);
end;

%% Plot data and autocorrelation function 

f=figure('name','Time series plot: Simulated data');
h=plot(Y(:,1),'k','linewidth',2);
hold on
h1=plot(Y(:,2),'color',[0.5 0.5 0.5],'linewidth',2);
set(gca,'xlim',[0 T-1]);
hl=legend([h h1],{'y1','y2'});
set(hl,'box','off','location','northeast')

for i=1:nvary
    for j=1:nvary
        plotautocorrTheory(['varacf_v' int2str(i) '_v' int2str(j)],permute(R(i,j,:),[3 1 2]));
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotautocorrTheory(name,acf)

hor=numel(acf);

f=figure('name','Autocorrelation Theory');
h=bar(acf,'k');
set(gca,'xtick',1:hor);        
axis tight
set(gca,'ylim',[-1 1]);















