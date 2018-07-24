clearvars; clc;close all;
% Data Handling
y = xlsread('../data/data.xlsx','RWZ2010');
% Options
opt.vnames   = {'GDP level','Federal Funds Rate', 'Deflator Inflation'}; 
opt.epsnames = {'Policy','Aggregate Demand','Aggregate Supply'}; %note that structural shocks are not ordered in the same row as the corresponding variable
opt.IRFcumsum = [1 0 0];
opt.filename = 'RWZ';
opt.nsteps = 40;
opt.nlag = 4;
opt.doplot = 1;
opt.dosave = 1;
% Estimate reduced-form
nlag = opt.nlag;
const = 1;
VAR = VARReducedForm(y,nlag,const);
A1inv_big = inv(eye(size(VAR.Acomp,1))-VAR.Acomp); % from the companion
%J = [eye(VAR.nvar),zeros(VAR.nvar,VAR.nvar*(VAR.nlag-1))];
%LRMat = J*A1inv_big*J';
LRMat = A1inv_big(1:VAR.nvar,1:VAR.nvar);    % total impact matrix inv(eye(nvars)-A1hat-A2hat-...-Aphat)

% Structural identification
% Identification restrictions are set up in RWZSRLR_f
f = 'RWZSRLR_f';
StartValueMethod = 2; %0: Use identity matrix, 1: use square root, 2: use cholesky as starting value 
% Options for fsolve
TolX = 1e-4;            % Termination tolerance on the current point
TolFun = 1e-9;          % Termination tolerance on the function value
MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
MaxIter = 1000;         % Maximum numberof iterations allowed
OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
% Initital guess
if StartValueMethod == 0
    B0inv = eye(VAR.nvar); % Use identity matrix as starting value
elseif StartValueMethod == 1
    B0inv = VAR.SIGu^.5; % Use square root of vcov of reduced form as starting value
elseif StartValueMethod == 2
    B0inv = chol(VAR.SIGu,'lower'); % Use Cholesky decomposition of vcov of reduced form
end
% Call optimization routine fsolve
[B0inv,fval,exitflag,output]=fsolve(f,B0inv,options,VAR.SIGu,LRMat);

% Normalize sign of first column such that a monetary policy shock (first column) raises the interest rate (second row) (monetary tightening)
if sign(B0inv(2,1)) == -1
    B0inv(:,1)=-B0inv(:,1);
end
% Normalize sign of second column such that a positive aggregate demand shock (2nd column) does not lower real GNP (first row) and inflation (third row)
if sign(B0inv(1,2)) == -1 && sign(B0inv(3,2)) == -1
    B0inv(:,2)=-B0inv(:,2);
end
% Normalize sign of third column such that a positive aggregate supply shock (3rd column) does not lower real GNP (first row) and does not raise inflation (third row)
if sign(B0inv(1,3)) == -1 && sign(B0inv(3,3)) == 1
    B0inv(:,3)=-B0inv(:,3);
end

table(B0inv)

% check that B0inv solution is correct (result should be (close to a) K x K zero matrix)
B0inv*B0inv'-VAR.SIGu
% check that structural innovations are orthogonal to one another (result should be identity matrix for correlations)
E=inv(B0inv)*VAR.residuals(1:VAR.nvar,:); 
corrcoef(E(2,:),E(1,:))

% Compute and plot structural impulse response function
IRFpoint = IRFs(VAR.Acomp,B0inv,opt);