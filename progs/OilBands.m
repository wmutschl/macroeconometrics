% Load monthly data for 1973.2-2007.12 ordered as:
% 1. Growth rate of world oil production 
% 2. Global real activity (business cycle index based on dry cargo shipping rates)
% 3. Real price of oil 
% The data sources are described in the text.
clearvars; clc;close all;
% Data Handling
ENDO = xlsread('usoilmonthly.xlsx',1);
ENDO=detrend(ENDO,0);
% Estimate reduced-form
nlag = 24;
const = 1;
opt.IRFcumsum = [1 0 0];
opt.vnames = {'Oil Production','Real Activity', 'Real Price of Oil'}; 
opt.epsnames = {'Oil Supply','Aggregate Demand', 'Oil-specific Demand'}; 

ENDO = xlsread('../data/data.xlsx','USOil');
% Estimate reduced-form
nlag = 4;
const = 1;
opt.IRFcumsum = [1 1 1];
opt.vnames = {'Real Price of Oil','GDP Deflator', 'Real GDP'}; 
%opt.epsnames = {'Oil Supply','Aggregate Demand', 'Oil-specific Demand'}; 

{'Real Price of Oil','GDP Deflator', 'Real GDP'};
VAR = VARReducedForm(ENDO,nlag,const);
% Structural identification
B0inv = chol(VAR.SIGu,'lower');
% Normalize sign of B0inv such that diagonal elements are positive
if sum(diag(B0inv)<0) ~= 0
    x = diag(B0inv)<0;
    B0inv(:,find(x==1)) = -1*B0inv(:,find(x==1));
end
% Consider oil price disruption
% Normalize sign of first column
B0inv(:,1)=-B0inv(:,1);
table(B0inv)
    
% Compute structural impulse response function
opt.filename = 'USOilMonthly';
opt.nsteps = 15;
opt.nlag = VAR.nlag;
opt.doplot = 0;
opt.dosave = 0;
IRFpoint = IRFs(VAR.Acomp,B0inv,opt);
% Compute bootstrap standard errors of IRFs
randn('seed',1234);
IRFse = bootstd(VAR,opt);

if 1<0
    % VAR bootstrap
    randn('seed',1234);
    nrep=200;
    trmat=zeros(VAR.nvar,VAR.nvar,opt.nsteps+1,r);
    parfor_progress(nrep); % Initialize 
    parfor r=1:nrep
        ystar = BootstrapGDP(VAR);
        VARstar = VARReducedForm(ystar,nlag,const,0);
        B0invstar = chol(VARstar.SIGu,'lower');
        IRFstar = IRFs(VARstar.Acomp,B0invstar,opt);
        IRFrsestar = bootstd(VARstar,opt);
        trmat(:,:,:,r)=(IRFstar-IRFpoint)./IRFse; % Symmetric percentile-t
        parfor_progress; % Count 
    end
    parfor_progress(0); % Clean up
    % Calculate 95 percent symmetric percentile t interval endpoints
    CILO=IRFpoint-prctile(trmat,97.5,4).*IRFse;
    CIUP=IRFpoint-prctile(trmat,2.5,4).*IRFse;
    CILO(isnan(CILO))=0; CIUP(isnan(CIUP))=0;
else
    CILO=IRFpoint-1.96*IRFse; CIUP=IRFpoint+1.96*IRFse;
end


figure('Name','Inference');
countplots = 1;
x_axis = zeros(1,opt.nsteps+1);
for ishock = 1:VAR.nvar
    for ivar = 1:VAR.nvar
        subplot(VAR.nvar,VAR.nvar,countplots);        
        irfpoint = squeeze(IRFpoint(ivar,ishock,:));
        irfup    = squeeze(CIUP(ivar,ishock,:));
        irflo    = squeeze(CILO(ivar,ishock,:));
        plot(0:1:opt.nsteps,irfpoint,'b','LineWidth',2);        
        hold on;        
        plot(0:1:opt.nsteps, [irflo irfup] ,'--r');
        plot(0:1:opt.nsteps,x_axis,'k','LineWidth',2);
        grid;
        xlim([0 opt.nsteps]);
        title(opt.vnames{ivar})
        ylabel([opt.epsnames{ishock}, 'Shock'])
        countplots = countplots + 1;
    end
end

