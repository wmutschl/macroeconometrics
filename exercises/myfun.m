function chapter1_gdpAndNormality()
%% PURPOSE: Simple distribution plots
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
data_raw=xlsread('nwgdp.xlsx','q')';

dates=data_raw(:,1);
data=data_raw(:,2);

%% Figures
f=figure('name','Histogram of data');
h=histfit(data,10);
set(h(2),'color',[0 0 0])
set(h(1),'facecolor',[0.5 0.5 0.5])

f=figure('name','Normplot of data');
h=normplot(data);
set(h(1),'markeredgecolor',[0 0 0]);
set(h(2),'color',[0.5 0.5 0.5],'linewidth',2)
set(h(3),'color',[0.5 0.5 0.5],'linewidth',2)
grid off
set(get(gca,'title'),'string','')
set(get(gca,'xlabel'),'string','')
set(get(gca,'ylabel'),'string','')
set(gca,'yticklabel',num2str(str2double(get(gca,'yticklabel'))*100,'%6.0f'))

f=figure('name','Data');
h=plot(data,'k','linewidth',2);

