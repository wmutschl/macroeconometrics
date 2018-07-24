% TRIVAR.M
% Lutz Kilian
% University of Michigan
% April 1997

%clear all;

global h t

% Load monthly data for 1973.2-2007.12 ordered as:
% 1. Growth rate of world oil production 
% 2. Global real activity (index based on dry cargo shipping rates)
% 3. Real price of oil 
% The data sources are described in the text.

load data.txt; y=data; [t,q]=size(y); 
time=(1973+2/12:1/12:2007+12/12)';  % Time line
h=15;                               % Impulse response horizon
p=24;                               % VAR lag order
y=detrend(y,0);
[A,SIGMA,Uhat,V,X]=olsvarc(y,p);	% VAR with intercept	
SIGMA=SIGMA(1:q,1:q);
