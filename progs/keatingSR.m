clearvars; clc; close all;
% Data Handling
y = xlsread('../data/data.xlsx','Keating1992');
% Estimate reduced-form
nlag = 4;
const = 1;
VAR = VARReducedForm(y,nlag,const);
% Identification restrictions are set up in KeatingSR_f_SR
f = 'KeatingSR_fSR';
% Options for fsolve
TolX = 1e-4;            % Termination tolerance on the current point
TolFun = 1e-9;          % Termination tolerance on the function value
MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
MaxIter = 1000;         % Maximum numberof iterations allowed
OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
% Initital guess, note that candidates are of the structure [B0 diag(SIGeps)]
B0_diageps= [eye(VAR.nvar) ones(VAR.nvar,1)];
% Call optimization routine fsolve to minimize f
KeatingSR_fSR(B0_diageps,VAR.SIGu)'
[B0_diageps,fval,exitflag,output]=fsolve(f,B0_diageps,options,VAR.SIGu);
B0 = B0_diageps(:,1:VAR.nvar);
SIGeps = B0_diageps(:,VAR.nvar+1); % These are just the variances
impact = inv(B0)*diag(sqrt(SIGeps)); % Compute impact matrix used to plot IRFs, i.e. inv(B0)*sqrt(SigEps)
table(impact)
    
% Compute and plot structural impulse response function
opt.filename = 'Keating1992';
opt.nsteps = 12;
opt.nlag = VAR.nlag;
opt.doplot = 1;
opt.dosave = 1;
opt.IRFcumsum = [1 1 0 1];
opt.vnames   = {'Deflator','Real GNP', 'Federal Funds Rate', 'M1'};
opt.epsnames = {'AS','IS','MS','MD'};
IRFpoint = IRFs(VAR.Acomp,impact,opt);