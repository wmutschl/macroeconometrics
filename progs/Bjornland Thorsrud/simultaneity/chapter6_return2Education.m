function [results_tsls,results_ols,results_olsn]=chapter6_return2Education()
% PURPOSE: Do instrumental variable regression (TSLS) to find effects of
% eduction on wage
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

%% Load data
fid = fopen('./asciiqob.txt','rt');
A = textscan(fid, '%f %f %f %f %f', 'HeaderLines',0);
fclose(fid);

% log_weekly_wage education year_of_birth quarter_of_birth place_of_birth 
log_weekly_wage=A{1};
education=A{2};
year_of_birth=A{3};
quarter_of_birth=A{4};
place_of_birth=A{5};

j=1:3;
c=1:10;
cind=c+(30-1);

N=size(log_weekly_wage,1);

% create born in quarter dummy
Q=zeros(N,max(j));
for i=1:N
    if quarter_of_birth(i)==1
        Q(i,1)=1;
    elseif quarter_of_birth(i)==2
        Q(i,2)=1;
    elseif quarter_of_birth(i)==3
        Q(i,3)=1;
    end;
end;
% create year dummy
Y=zeros(N,max(c));
for i=1:N
    idx=year_of_birth(i)==cind;        
    Y(i,idx)=1;
end;

% Create interacting variable variable
QY=zeros(N,max(c)*max(j));
for ii=1:N
    cnt=1;
    for jj=1:max(c)
        for kk=1:max(j)
            QY(ii,cnt)=Y(ii,jj)*Q(ii,kk);
            cnt=cnt+1;
        end;
    end;    
end;

% Constant
C=ones(N,1);

%% TSLS

% Call TSLS function (implements TSLS differently than in the book, but
% easier to compute std etc.) 
results_tsls=tslsEstimation(log_weekly_wage,education,[C Y(:,1:end-1)],[C Y(:,1:end-1) QY]);
% More transparent formula:
% % First stage
%     x=[C Y(:,1:end-1) QY];
%     beta1=x\education;
%     education_hat=x*beta1;
% % Second stage
%     beta_tsls=[education_hat C Y(:,1:end-1)]\log_weekly_wage;

%% OLS
results_ols=olsEstimation(log_weekly_wage,[education C Y(:,1:end-1)]);

%% OLS naive
results_olsn=olsEstimation(log_weekly_wage,[education C]);






