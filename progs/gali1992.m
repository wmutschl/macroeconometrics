clearvars; clc;close all;
% Data Handling
y = xlsread('gali1992.xlsx');
y = detrend(y,'constant');
% Options
opt.vnames   = {'dgnp','di', 'i-dp', 'dm-dp'}; 
opt.epsnames = {'AS','MS','MD','IS'};
opt.IRFcumsum = [1 1 0 0];
opt.filename = 'gali1992';
opt.nsteps = 40;
opt.nlag = 4;
opt.doplot = 1;
opt.dosave = 1;
% Estimate reduced-form
nlag = opt.nlag;
const = 1;
VAR = VARReducedForm(y,nlag,const,0);
VAR.SIGu
VAR.nobse

A1inv_big = inv(eye(size(VAR.Acomp,1))-VAR.Acomp); % from the companion
J = [eye(VAR.nvar),zeros(VAR.nvar,VAR.nvar*(VAR.nlag-1))];
LRMat = J*A1inv_big*J';

% Structural identification
% Identification restrictions are set up in f
f = 'gali1992_f';
StartValueMethod = 2; %0: Use identity matrix, 1: use square root, 2: use cholesky as starting value 
% Options for fsolve
TolX = 1e-4;            % Termination tolerance on the current point
TolFun = 1e-9;          % Termination tolerance on the function value
MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
MaxIter = 1000;         % Maximum numberof iterations allowed
OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
%OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve

options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
% Initital guess
if StartValueMethod == 0
    B0inv = eye(VAR.nvar); % Use identity matrix as starting value
elseif StartValueMethod == 1
    B0inv = VAR.SIGu^.5; % Use square root of vcov of reduced form as starting value
elseif StartValueMethod == 2
    B0inv = chol(VAR.SIGu,'lower'); % Use Cholesky decomposition of vcov of reduced form
end
[B0inv,fval,exitflag,output]=fsolve(f,B0inv,options,VAR.SIGu,LRMat);
B0inv
LRMat*B0inv

% % From EVIEWS Ae=Bu with E(uu')=I
A = [1 -0.655408 0.293322 -0.086442;
    0.239160 1 0.301674 -0.301674;
    -1.557400 4.985212 1 1.213470;
    4.654613 7.582059 -3.393284 1];
B = [1.070711 0 0 0;
    0 1.699136 0 0;
    0 0 6.515471 0;
    0 0 0 8.202743];
inv(A)*B
B0inv
% check that B0inv solution is correct (result should be (close to a) K x K zero matrix)
B0inv*B0inv'-VAR.SIGu
% check that structural innovations are orthogonal to one another (result should be identity matrix for correlations)
E=inv(B0inv)*VAR.residuals(1:VAR.nvar,:); 
corrcoef(E(2,:),E(1,:))

% Compute and plot structural impulse response function
IRFpoint = IRFs(VAR.Acomp,B0inv,opt);