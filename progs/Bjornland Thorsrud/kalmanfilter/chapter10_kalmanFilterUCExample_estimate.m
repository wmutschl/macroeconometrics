function result=chapter10_kalmanFilterUCExample_estimate()
% PURPOSE: Estimate UC model using Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The code is not written for computational efficiency or
% elegance. The book: 
%
% "Applied Time Series for Macroeconomics"
% Gyldendal Akademisk 2015
% by Hilde C. Bjørnland and Leif A. Thorsrud 
%
% provides details. Please refer to the book if the code(s) are used for 
% research of commercial purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings

% Which observation to start evaluating ML at (in maximization rutine)
start=5;

%% Load data

[data,~]=xlsread('.\usgdp.xlsx','q');
dates=data(1,:)';
y=data(2,:)';
t=size(y,1);
 
%% Kalman Filter ML

fun=@(x)KFoptim(x,y,start);       

% x0 initial opt values: ma ar c sigma
x0=zeros(5,1);
x0(end-1)=1;
x0(end)=1;

% Optimization settings
optOptions=optimset('TolX',1e-10,'TolFun',1e-10,'Maxiter',25000,'Display','iter','MaxFunEvals',25000);    
% Do optimization
[x1_unc,fval_unc,exitflag_unc,optOutput_unc,grad_unc,hessian_unc]=fminunc(fun,x0,optOptions);
% Get final estimated states and log likelihood evaluated over whole sample
[state_unc, varCov_unc, loglik_unc, ap_unc, Pp_unc,smoothState_unc,smoothVarCov_unc]=KFfinal(x1_unc,y);

%% Transform output
z1=x1_unc(1)/(1+abs(x1_unc(1)));
z2=x1_unc(2)/(1+abs(x1_unc(2)));
x1(1)=z1+z2;
x1(2)=-1*z1*z2;
x1(3)=x1_unc(3);
x1(4)=exp(x1_unc(4));
x1(5)=exp(x1_unc(5));

%%
result.b=x1(1:end-2);
result.sigma2=x1(end-1:end);
result.loglik=sum(loglik_unc);
result.trend=permute(smoothState_unc(1,1,:),[3 1 2]);
result.cycle=permute(smoothState_unc(2,1,:),[3 1 2]);
result.y=y;
result.dates=dates;
result.start=start;
result.nobs=size(y,1);

%% HP as UC
[state, varCov, loglik, ap, Pp,smoothState,smoothVarCov]=HPasUC(y);

result_hpasuc.loglik=sum(loglik);
result_hpasuc.trend=permute(smoothState(1,1,:),[3 1 2]);
result_hpasuc.cycle=y-result_hpasuc.trend;
%% Output figures and tables

ps=printSettings();

f=figure('name','Trend and actual');
h=plot(y,'k','linewidth',ps.lineWidth);
hold on
h1=plot(result.trend,'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','-');
axis tight
lh=legend([h h1],{'Actual','UC trend'});
set(lh,'box','off','location','northwest');
printSettings(gca,f);

f=figure('name','Cycle');
h=plot(result.cycle,'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth);
hold on
h1=plot(zeros(result.nobs,1),'k','linewidth',ps.lineWidth);
h2=plot(result_hpasuc.cycle,'color',[0 0 0],'linewidth',ps.lineWidth);
axis tight
lh=legend([h h2],{'UC cycle','HP as UC'});
set(lh,'box','off','location','northwest');
printSettings(gca,f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LL=KFoptim(x,y,start)
% Purpose: Runs the kalman filter with input variables. Returns the sum of
% the log likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Observation loadings
Z=[1 1 0];
% Observation constant;
d=0;
% Observatioin error (=0, no measurement error)
H=0;

z1=x(1)/(1+abs(x(1)));
z2=x(2)/(1+abs(x(2)));

phi1=z1+z2;
phi2=-1*z1*z2;

% Transition matrix
T=[1 0 0;
    0 phi1 phi2;
    0 1 0];
% Transition error (take exp() to ensure positive)
Q=zeros(2);
Q(1,1)=exp(x(4));
Q(2,2)=exp(x(5));
% Transition constant (zeros)
c=zeros(3,1);
c(1)=x(3);
% Matrix to map Q to states
R=zeros(3,2);
R(1:2,1:2)=eye(2);
% Initial values
a0=zeros(3,1);
a0(1)=y(1);
P0=eye(3)*100000;

[~,~,loglik,~,~]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
LL=-sum(loglik(start:end));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state, varCov, loglik, ap, Pp,smoothState,smoothVarCov]=KFfinal(x,y)
% PURPOSE: Run the Kalman filter with optimized values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=size(y,1);

% Observation loadings
Z=[1 1 0];
% Observation constant;
d=0;
% Observatioin error (=0, no measurement error)
H=0;

z1=x(1)/(1+abs(x(1)));
z2=x(2)/(1+abs(x(2)));

phi1=z1+z2;
phi2=-1*z1*z2;

% Transition matrix
T=[1 0 0;
    0 phi1 phi2;
    0 1 0];
% Transition error (take exp() to ensure positive)
Q=zeros(2);
Q(1,1)=exp(x(4));
Q(2,2)=exp(x(5));
% Transition constant (zeros)
c=zeros(3,1);
c(1)=x(3);
% Matrix to map Q to states
R=zeros(3,2);
R(1:2,1:2)=eye(2);
% Initial values
a0=zeros(3,1);
P0=eye(3)*100;

[state,varCov,loglik,ap,Pp]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
[smoothState,smoothVarCov]=KalmanSmoother(y',state,varCov,ap,Pp,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state, varCov, loglik, ap, Pp,smoothState,smoothVarCov]=HPasUC(y)
% PURPOSE: Run the Kalman filter with optimized values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=size(y,1);

% Observation loadings
Z=[1 0];
% Observation constant;
d=0;
% Observatioin error (=0, no measurement error)
H=1600;

% Transition matrix
T=[1 1;
   0 1];
% Transition error (take exp() to ensure positive)
Q=1; 
% Transition constant (zeros)
c=zeros(2,1);
% Matrix to map Q to states
R=zeros(2,1);
R(2,1)=1;
% Initial values
a0=zeros(2,1);
a0(1)=y(1);
P0=eye(2)*10000;

[state,varCov,~,ap,Pp]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
[smoothState,smoothVarCov]=KalmanSmoother(y',state,varCov,ap,Pp,T);

a0=smoothState(:,1,1);
P0=smoothVarCov(:,:,1);
[~,~,loglik,~,~]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
