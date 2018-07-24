clear; close all;
y = xlsread('../data/data.xlsx','ThreeVariableVAR');
varnames = {'Real GNP Growth' 'Federal Funds Rate' 'GNP Deflator Inflation'};
dates = (1954.75:0.25:2007.75)';
for j=1:3
    subplot(3,1,j); 
    plot(dates,y(:,j));    
    title(varnames{j});
end
nlag=4;
const=1;
VAR4 = VARReducedForm(y,nlag,const);
VAR4.maxEig
for j=1:VAR4.nvar
    disp(transpose(VAR4.(sprintf('eq%d',j)).tprob) <= 0.05);
end