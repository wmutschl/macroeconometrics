% =======================================================================
% Replicates the Blanchard and Quah (1989) model using long-run restrictions 
% to identify the structural shocks either via a Cholesky decomposition or
% using a numerical solver for expository purposes
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
% =======================================================================

clearvars; clc;close all;

%% Settings and options
identifmethod = 1; % 0: use Cholesky decomposition, 1: use fsolve to numerically find the impact matrix

const = 1;         % 0: no constant, 1: constant, 2: constant and linear trend
nlag = 8;          % Number of lags

opt.nlag = nlag;                             % Number of lags
opt.filename = 'figures/BlanchardQuah1989';          % Filename to save results
opt.nsteps = 40;                             % Horizon of IRFs
opt.doplot = 1;                              % 1: Plot IRF functions
opt.dosave = 1;                              % 1: Save IRF functions to file
opt.IRFcumsum = [1 0];                       % Cumulate (1) or not (0) IRFs for each variable
opt.vnames   = {'GDP level','Unemployment'}; % Variable names in IRF plots
opt.epsnames = {'Supply','Demand'};          % Shock names in IRF plots

if identifmethod == 1
    % File where identification restrictions are set up
    f = 'BlanchardQuah1989_f';
    StartValueMethod = 1;   % 0: Use identity matrix, 1: use square root, 2: use cholesky as starting value 
    % Options for fsolve
    TolX = 1e-4;            % Termination tolerance on the current point
    TolFun = 1e-9;          % Termination tolerance on the function value
    MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
    MaxIter = 1000;         % Maximum numberof iterations allowed
    OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
    options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm); % Make option set
end

%% Data Handling
ENDO = xlsread('../data/data.xlsx','BlanchardQuah1989'); % Load data

%% Reduced-form Estimation
VAR = VARReducedForm(ENDO,nlag,const); % Estimate reduced form
A1inv_big = inv(eye(size(VAR.Acomp,1))-VAR.Acomp); % Long-run matrix using companion form
J = [eye(VAR.nvar),zeros(VAR.nvar,VAR.nvar*(VAR.nlag-1))];
LRMat = J*A1inv_big*J';  % total impact matrix inv(eye(nvars)-A1hat-A2hat-...-Aphat)

%% Structural Estimation
if identifmethod == 0
    THETA = chol(LRMat*VAR.SIGu*LRMat','lower');
    B0inv = inv(LRMat)*THETA;
elseif identifmethod == 1
    if StartValueMethod == 0
        B0inv = eye(VAR.nvar);          % Use identity matrix as starting value
    elseif StartValueMethod == 1
        B0inv = VAR.SIGu^.5;            % Use square root of vcov of reduced form as starting value
    elseif StartValueMethod == 2
        B0inv = chol(VAR.SIGu,'lower'); % Use Cholesky decomposition of vcov of reduced form
    end
    % Call optimization routine fsolve to minimize f
    [B0inv,fval,exitflag,output]=fsolve(f,B0inv,options,VAR.SIGu,LRMat);
end

% Normalization rule on impact matrix: Diagonal elements of B0inv are supposed to be positive
if sum(diag(B0inv)<0) ~= 0 
    x = diag(B0inv)<0;
    B0inv(:,find(x==1)) = -1*B0inv(:,find(x==1)); 
end
table(B0inv)

%% Compute and plot structural impulse response function
IRFpoint = IRFs(VAR.Acomp,B0inv,opt);