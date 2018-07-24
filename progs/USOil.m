clearvars; clc;close all;
identifmethod = 1; % 0: for Cholesky 1: for fsolve
% Data Handling
y = xlsread('../data/data.xlsx','USOil');
% Estimate reduced-form
nlag = 4;
const = 1;
VAR = VARReducedForm(y,nlag,const);
% Structural identification
if identifmethod == 0
    B0inv = chol(VAR.SIGu,'lower');
elseif identifmethod == 1
    % Identification restrictions are put into USOil_fSR
    f = 'USOil_fSR';
    StartValueMethod = 1; %0: Use identity matrix, 1: use square root, 2: use cholesky as starting value 
    % Options for fsolve
    TolX = 1e-4;            % Termination tolerance on the current point
    TolFun = 1e-9;          % Termination tolerance on the function value
    MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
    MaxIter = 1000;         % Maximum numberof iterations allowed
    OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
    options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
    if StartValueMethod == 0
        B0inv = eye(VAR.nvar); % Use identity matrix as starting value
    elseif StartValueMethod == 1
        B0inv = VAR.SIGu^.5; % Use square root of vcov of reduced form as starting value
    elseif StartValueMethod == 2
        B0inv = chol(VAR.SIGu,'lower'); % Use Cholesky decomposition of vcov of reduced form
    end
    [B0inv,fval,exitflag,output]=fsolve(f,B0inv,options,VAR.SIGu);
end

% Normalize sign of B0inv such that diagonal elements are positive
if sum(diag(B0inv)<0) ~= 0
    x = diag(B0inv)<0;
    B0inv(:,find(x==1)) = -1*B0inv(:,find(x==1));
end
table(B0inv)
    
% Compute and plot structural impulse response function
opt.filename = 'USOil';
opt.nsteps = 12;
opt.nlag = VAR.nlag;
opt.doplot = 1;
opt.dosave = 1;
opt.IRFcumsum = [1 1 1];
opt.vnames = {'Real Price of Oil','GDP Deflator', 'Real GDP'}; 
IRFpoint = IRFs(VAR.Acomp,B0inv,opt);