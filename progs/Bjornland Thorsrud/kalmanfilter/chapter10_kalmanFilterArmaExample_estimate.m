function result=chapter10_kalmanFilterArmaExample_estimate()
% PURPOSE: Estimate ARMA model using Kalman Filter
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
% ARMA AR(p) lags
p=1;
% ARMA MA(q) lags
q=1;

% Optimization settings
optOptions=optimset('TolX',1e-10,'TolFun',1e-10,'Maxiter',25000,'Display','iter','MaxFunEvals',25000);    

%% Load data
[data,~]=xlsread('.\ukinf.xlsx','q');
dates=data(1,:)';
y=data(2,:)';
t=size(y,1);
 
%% Kalman Filter ML

fun=@(x)KFoptim(x,y,p,q,start);       
% x0 initial opt values: ma ar c sigma
x0=zeros(p+q+2,1);
x0(end-1)=mean(y);
x0(end)=1;
% Do optimization
[x1_unc,fval_unc,exitflag_unc,optOutput_unc,grad_unc,hessian_unc]=fminunc(fun,x0,optOptions);
% Get final estimated states and log likelihood evaluated over whole sample
[state_unc, varCov_unc, loglik_unc, ap_unc, Pp_unc,smoothState_unc,smoothVarCov_unc]=KFfinal(x1_unc,y,p,q);

%% Compute standard errors and p-values, R2 etc. unc
r=max(p,q+1);

bstd_unc=sqrt(diag(inv(hessian_unc)));

T=zeros(r,r);
T(1,1:p)=x1_unc(q+1:q+p);
T(2:end,1:end-1)=eye(r-1);
resid=(permute(smoothState_unc(:,1,2:end),[1 3 2])-T*permute(smoothState_unc(:,1,1:end-1),[1 3 2]))';
resid=resid(r:end,1);
yhat=y(r+1:end)-resid;

nobs=size(resid,1);
nvarx=p+q+1;    
betat=x1_unc(:)./bstd_unc(:);
pval=2*tcdf(-abs(betat),nobs-nvarx);                
m0=eye(nobs)-1/nobs.*ones(nobs);   
ee=diag(resid'*resid);                        
R2=1-ee/(y(r+1:end)'*m0*y(r+1:end));                        

result.b=x1_unc(1:end-1);
% remeber to transform sigma from logs
result.sigma2=exp(x1_unc(end));
result.loglik=sum(loglik_unc(r+1:end));
result.bstd=bstd_unc(1:end-1);
result.R2=R2;
result.pval=pval(1:end-1);
result.yhat=yhat;
result.resid=resid;

result.start=start;
result.p=p;
result.q=q;
result.nobs=nobs;

%% AR(1) for comparison

p=1;
q=0;

fun=@(x)KFoptim(x,y,p,q,start);       
% x0 initial opt values: ma ar c sigma
x0=zeros(p+q+2,1);
x0(end-1)=mean(y);
x0(end)=1;
% Do optimization
[x1_unc,fval_unc,exitflag_unc,optOutput_unc,grad_unc,hessian_unc]=fminunc(fun,x0,optOptions);
% Get final estimated states and log likelihood evaluated over whole sample
[state_unc, varCov_unc, loglik_unc, ap_unc, Pp_unc,smoothState_unc,smoothVarCov_unc]=KFfinal(x1_unc,y,p,q);

%% Compute standard errors and p-values, R2 etc. unc
r=max(p,q+1);

bstd_unc=sqrt(diag(inv(hessian_unc)));

T=zeros(r,r);
T(1,1:p)=x1_unc(q+1:q+p);
T(2:end,1:end-1)=eye(r-1);
resid=(permute(smoothState_unc(:,1,2:end),[1 3 2])-T*permute(smoothState_unc(:,1,1:end-1),[1 3 2]))';
resid=resid(r:end,1);
yhat=y(r+1:end)-resid;

nobs=size(resid,1);
nvarx=p+q+1;    
betat=x1_unc(:)./bstd_unc(:);
pval=2*tcdf(-abs(betat),nobs-nvarx);                
m0=eye(nobs)-1/nobs.*ones(nobs);   
ee=diag(resid'*resid);                        
R2=1-ee/(y(r+1:end)'*m0*y(r+1:end));                        

result_ar.b=x1_unc(1:end-1);
% remeber to transform sigma from logs
result_ar.sigma2=exp(x1_unc(end));
result_ar.loglik=sum(loglik_unc(r+1:end));
result_ar.bstd=bstd_unc(1:end-1);
result_ar.R2=R2;
result_ar.pval=pval(1:end-1);
result_ar.yhat=yhat;
result_ar.resid=resid;

result_ar.start=start;
result_ar.p=p;
result_ar.q=q;
result_ar.nobs=nobs;

%% Output figures and tables

ps=printSettings();

f=figure('name','Fitted and actual');
h=plot(y(r+1:end).*100,'k','linewidth',ps.lineWidth);
hold on
h1=plot(yhat.*100,'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','--');
axis tight
lh=legend([h h1],{'Inflation (QoQ)','Fitted values'});
set(lh,'box','off','location','northeast');
printSettings(gca,f);

f=figure('name','Residuals');
h=plot(resid.*100,'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth);
hold on
h1=plot(zeros(nobs,1),'k','linewidth',ps.lineWidth);
axis tight
printSettings(gca,f);

fid=1;
fprintf(fid,'\\begin{table}[t!]\\footnotesize \n');
fprintf(fid,'\\caption{Estimation results: AR(1) and ARMA(1,1) with UK inflation (QoQ)} \n');
fprintf(fid,'\\label{chap:kalmanfilter:tab:estimationResultsARMAukinf} \n');
fprintf(fid,'\\begin{center} \n');
fprintf(fid,'\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llccc} \\\\ \\toprule \n');
fprintf(fid,' &  & \\multicolumn{3}{l}{\\textbf{Parameters}}\\\\ \\cmidrule(r){3-5} \n');
fprintf(fid,'\\textbf{Model} & &\\textbf{$\\mu$} &\\textbf{$\\phi$}&\\textbf{$\\theta$} \\\\ \\midrule \n');
fprintf(fid,'\\textbf{ARMA(1,1)} & Estimate &  %6.4f &  %6.4f &  %6.4f \\\\ \n',result.b(3),result.b(2),result.b(1));
fprintf(fid,'          & p-value  &  %6.4f &  %6.4f &  %6.4f \\\\ \n',result.pval(3),result.pval(2),result.pval(1));
fprintf(fid,'          & $R^{2}$  & %6.2f & & \\\\ \n',result.R2);
fprintf(fid,'          & Log likelihood  & %6.2f & & \\\\ [0.25cm] \n',result.loglik);
fprintf(fid,'\\textbf{AR(1)} & Estimate &  %6.4f &  %6.4f &   \\\\ \n',result_ar.b(2),result_ar.b(1));
fprintf(fid,'      & p-value  &  %6.4f &  %6.4f &   \\\\ \n',result_ar.pval(2),result_ar.pval(1));
fprintf(fid,'          & $R^{2}$  & %6.2f & & \\\\ \n',result_ar.R2);
fprintf(fid,'          & Log likelihood  & %6.2f & & \\\\ [0.25cm] \n',result_ar.loglik);
fprintf(fid,'\\textbf{Number of observations} & %6.0f & & & \\\\ \n',nobs);
fprintf(fid,'\\textbf{Sample}   & \\multicolumn{2}{l}{1988:Q3 - 2013:Q3}  & & \\\\ \\bottomrule \n');
fprintf(fid,'\\end{tabular*}\n');
fprintf(fid,'\\end{center} \n');
fprintf(fid,'\\vspace{0.0cm}\n');
fprintf(fid,'\\begin{minipage}[!h]{\\textwidth}\n');
fprintf(fid,'\\textit{Note: See the text for details.}\n');
fprintf(fid,'\\end{minipage}\n');
fprintf(fid,'\\end{table} \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LL=KFoptim(x,y,p,q,start)
% Purpose: Runs the kalman filter with input variables. Returns the sum of
% the log likelihood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Totol number of ar and ma parameters
r=max(p,q+1);

% Observation loadings
Z=zeros(1,r);
Z(1:q+1)=[1 x(1:q)'];
% Observation constant;
d=x(end-1);
% Observatioin error (=0, no measurement error)
H=0;
% Transition matrix
T=zeros(r,r);
T(1,1:p)=x(q+1:q+p);
T(2:end,1:end-1)=eye(r-1);
% Transition error (take exp() to ensure positive)
Q=exp(x(end));
% Transition constant (zeros)
c=zeros(r,1);
% Matrix to map Q to states
R=zeros(r,r);
R(1,1)=1;
% Initial values
a0=zeros(r,1);
P0=eye(r)*100;

[~,~,loglik,~,~]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
LL=-sum(loglik(start:end));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [state, varCov, loglik, ap, Pp,smoothState,smoothVarCov]=KFfinal(x,y,p,q)
% PURPOSE: Run the Kalman filter with optimized values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Totol number of ar and ma parameters
r=max(p,q+1);
t=size(y,1);

% Observation loadings
Z=zeros(1,r);
Z(1:q+1)=[1 x(1:q)'];
% Observation constant;
d=x(end-1);
% Observatioin error (=0, no measurement error)
H=0;
% Transition matrix
T=zeros(r,r);
T(1,1:p)=x(q+1:q+p);
T(2:end,1:end-1)=eye(r-1);
% Transition error (take exp() to ensure positive)
Q=exp(x(end));
% Transition constant (zeros)
c=zeros(r,1);
% Matrix to map Q to states
R=zeros(r,r);
R(1,1)=1;
% Initial values
a0=zeros(r,1);
P0=eye(r)*100;

[state,varCov,loglik,ap,Pp]=KalmanFilter(y',Z,d,H,T,c,R,Q,a0,P0);
[smoothState,smoothVarCov]=KalmanSmoother(y',state,varCov,ap,Pp,T);