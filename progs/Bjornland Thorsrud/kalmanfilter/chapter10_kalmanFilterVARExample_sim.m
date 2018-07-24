function chapter10_kalmanFilterVARExample_sim()
% PURPOSE: Simulate VAR data, filter VAR data when data is missing
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

%% Simulate VAR data
[y,phi,sigma2,mu]=simulateVAR();

[t,nobs]=size(y);

%% Kalman Filter

% make KF input
Z=eye(nobs);
d=zeros(nobs,1);
H=zeros(nobs);
T=phi;
Q=sigma2;
c=mu;
R=eye(nobs);

a0=zeros(nobs,1);
P0=zeros(nobs);

%% Kalman Filter, now with missing values 1

yy=y;
yy([(10:30) (80:90)],1)=nan;

[state_1,varCov_1,loglik,ap_1,Pp_1]=KalmanFilter(yy',Z,d,H,T,c,R,Q,a0,P0);
[smoothState_1,smoothVarCov_1]=KalmanSmoother(yy',state_1,varCov_1,ap_1,Pp_1,T);

%% Kalman Filter, now with missing values 2

yy=y;
yy(:,1)=nan;

[state_2,varCov_2,loglik_2,ap_2,Pp_2]=KalmanFilter(yy',Z,d,H,T,c,R,Q,a0,P0);
[smoothState_2,smoothVarCov_2]=KalmanSmoother(yy',state_2,varCov_2,ap_2,Pp_2,T);

%% Kalman Filter, now with missing values 3

yy=y;
yy(t-10:t,1)=nan;

[state_3,varCov_3,loglik_3,ap_3,Pp_3]=KalmanFilter(yy',Z,d,H,T,c,R,Q,a0,P0);
[smoothState_3,smoothVarCov_3]=KalmanSmoother(yy',state_3,varCov_3,ap_3,Pp_3,T);

%% Kalman Filter, now with missing values 4

yy=y;
yy(t-10:t,1)=nan;
yy(t-10:t,2)=nan;

[state_4,varCov_4,loglik_4,ap_4,Pp_4]=KalmanFilter(yy',Z,d,H,T,c,R,Q,a0,P0);
[smoothState_4,smoothVarCov_4]=KalmanSmoother(yy',state_4,varCov_4,ap_4,Pp_4,T);

%% Make figures for output

ps=printSettings();

f=figure('name','KF - missing values');
h=plot(y(:,1),'color',[0 0 0],'linewidth',ps.lineWidth,'linestyle','-','marker','*');
hold on
h1=plot(permute(state_1(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','--');
h2=plot(permute(smoothState_1(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','-');
ylim=get(gca,'ylim');
line([t-11 t-11],[ylim(1) ylim(2)],'linewidth',2,'color',[0 0 0]);
set(gca,'xlim',[0 t-1]);
hl=legend([h h1 h2],{'Actual','Filtered','Smoothed'});
set(hl,'box','off','location','northwest')
printSettings(gca,f);

f=figure('name','KF - missing vector');
h=plot(y(:,1),'color',[0 0 0],'linewidth',ps.lineWidth,'linestyle','-','marker','*');
hold on
h1=plot(permute(state_2(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','--');
h2=plot(permute(smoothState_2(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','-');
ylim=get(gca,'ylim');
line([t-11 t-11],[ylim(1) ylim(2)],'linewidth',2,'color',[0 0 0]);
set(gca,'xlim',[0 t-1]);
hl=legend([h h1 h2],{'Actual','Filtered','Smoothed'});
set(hl,'box','off','location','northwest')
printSettings(gca,f);

f=figure('name','KF - conditional forecast');
h=plot([y(1:end-11,1);nan(10,1)],'color',[0 0 0],'linewidth',ps.lineWidth,'linestyle','-','marker','*');
hold on
h1=plot(permute(state_3(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','--');
h2=plot(permute(smoothState_3(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','-');
ylim=get(gca,'ylim');
line([t-11 t-11],[ylim(1) ylim(2)],'linewidth',2,'color',[0 0 0]);
set(gca,'xlim',[0 t-1]);
hl=legend([h h1 h2],{'Actual','Filtered','Smoothed'});
set(hl,'box','off','location','northwest')
printSettings(gca,f);

f=figure('name','KF - forecast');
h=plot([y(1:end-11,1);nan(10,1)],'color',[0 0 0],'linewidth',ps.lineWidth,'linestyle','-','marker','*');
hold on
h1=plot(permute(state_4(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','--');
h2=plot(permute(smoothState_4(1,1,:),[3 2 1]),'color',[0.5 0.5 0.5],'linewidth',ps.lineWidth,'linestyle','-');
ylim=get(gca,'ylim');
line([t-11 t-11],[ylim(1) ylim(2)],'linewidth',2,'color',[0 0 0]);
set(gca,'xlim',[0 t-1]);
hl=legend([h h1 h2],{'Actual','Filtered','Smoothed'});
set(hl,'box','off','location','northwest')
printSettings(gca,f);

f=figure('name','KF - state 2');
h=plot(y(:,2),'color',[0 0 0],'linewidth',ps.lineWidth,'linestyle','-','marker','*');
ylim=get(gca,'ylim');
line([t-11 t-11],[ylim(1) ylim(2)],'linewidth',2,'color',[0 0 0]);
set(gca,'xlim',[0 t-1]);
hl=legend(h,{'Observable'});
set(hl,'box','off','location','northwest')
printSettings(gca,f);

