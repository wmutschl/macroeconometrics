function result=chapter1_opennessAndGrowthOLS()
%% PURPOSE: Estimate relationship between GDP growth and openness 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The code is not written for computational efficiency or
% elegance. The book: 
%
% "Applied Time Series for Macroeconomics"
% Gyldendal Akademisk 2014
% by Hilde C. Bj�rnland and Leif A. Thorsrud 
%
% provides details. Please refer to the book if the code(s) are used for 
% research of commercial purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
[data,~]=xlsread('gdpAndOpenness.xlsx','DATA');
gdp=data(:,1);
openc=data(:,2);

%% Plot the data
f=figure('name','Scatter of data');
h=scatter(openc,gdp,'fill');
h1=lsline;
set(h,'CData',[0 0 0])
set(h1,'LineWidth',2,'Color',[0 0 0])

%% Do OLS 
%(Already a OLS line in the plot, but let us compute the coeffs., 
% t-stat etc.. Note, in the code below this is done using more or less the 
% same notation etc. as in Chapter 1 of the book. This is not efficient
% etc. but done for clarity.)

y=gdp;
x=openc;
% Number of observations
n=numel(y);
% Do OLS
bary=mean(y);
barx=mean(x);
hatbeta1=(x-barx)'*(y-bary)/sum((x-barx).^2);
hatbeta0=bary-hatbeta1*barx;
% Note: In Matlab OLS is easier
% hatbeta=[ones(n,1) x]\y
haty=hatbeta0+hatbeta1.*x;
hatu=y-haty; 

EES=sum((haty-bary).^2);
TSS=sum((y-bary).^2);

% Measures of fit
R2=EES/TSS;
SER=sqrt(sum(hatu.^2)/(n-2));
s2=SER^2;
xxinv=diag(inv([ones(n,1) x]'*[ones(n,1) x]));
SE=sqrt(s2.*xxinv);

% t-statistics
tstatb0=hatbeta0/SE(1);
tstatb1=hatbeta1/SE(2);
z2=abs(norminv(0.05/2,0,1));

% p-values
pvalb1=2*tcdf(-abs(tstatb1),n-2);                
pvalb0=2*tcdf(-abs(tstatb0),n-2);                


%% Put (some stuff) into result struct for output

result.pvalb1=pvalb1;
result.pvalb0=pvalb0;

result.tstatb0=tstatb0;
result.tstatb1=tstatb1;

result.R2=R2;

result.hatbeta0=hatbeta0;
result.hatbeta1=hatbeta1;










