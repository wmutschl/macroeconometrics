function chapter5_spectralDecompostitions() 
% PURPOSE: Compute and illustrate spectral decompositions
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

%% Generate data 2
T=100;
t=(1:1:T)';
w=[6/T 30/T 40/T];
x1=2*cos(2*pi*t*w(1))+3*sin(2*pi*t*w(1)); % frequency 6 cycles in T points 
x2=4*cos(2*pi*t*w(2))+5*sin(2*pi*t*w(2)); % frequency 10 cycles in T points 
x3=6*cos(2*pi*t*w(3))+7*sin(2*pi*t*w(3)); % frequency 40 cycles in T points 
y=x1+x2+x3;

%% Plot data

f=figure('name','Time series plot x1');
h=plot(x1,'k','linewidth',2);

f=figure('name','Time series plot x2');
h=plot(x2,'k','linewidth',2);

f=figure('name','Time series plot x3');
h=plot(x3,'k','linewidth',2);

f=figure('name','Time series plot y');
h=plot(y,'k','linewidth',2);

%% Do periodogram with plot
[Y,F]=spectralDecomp(y);
f=figure('name','Periodogram');
h=plot(F,Y,'k','linewidth',2,'lineStyle','--','marker','*');

[Y,F]=spectralDecomp(x1);
f=figure('name','Periodogram');
h=plot(F,Y,'k','linewidth',2,'lineStyle','--','marker','*');

[Y,F]=spectralDecomp(x2);
f=figure('name','Periodogram');
h=plot(F,Y,'k','linewidth',2,'lineStyle','--','marker','*');

[Y,F]=spectralDecomp(x3);
f=figure('name','Periodogram');
h=plot(F,Y,'k','linewidth',2,'lineStyle','--','marker','*');

%% Take away stuff between soem given frequencies using BP filter
pl=16; 
pu=20; 
[ybp,~]=bpass(y,pl,pu);
[Y_bp,F_bp]=spectralDecomp(ybp);

f=figure('name','Periodogram BP filtered');
h=plot(F_bp,Y_bp,'k','linewidth',2,'lineStyle','--','marker','*');

f=figure('name','Time series plot x1 and ybp');
h=plot(x1,'k','linewidth',2);
hold on
h1=plot(ybp,'color',[0.5 0.5 0.5],'linewidth',2);
set(gca,'ylim',[-4 5]);
hl=legend([h h1],{'Actual','Filtered'});
set(hl,'box','off','location','northwest')

f=figure('name','Time series plot x1 and ybp');
h=plot(y,'color',[0.5 0.5 0.5],'linewidth',2);
hold on
h1=plot(ybp,'k','linewidth',2);
hl=legend([h h1],{'Actual','Filtered'});
set(hl,'box','off','location','northwest')



