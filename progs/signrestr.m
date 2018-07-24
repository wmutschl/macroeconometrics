% =======================================================================
% Estimating a monetary VAR with 3 variables: real GDP, CPI, federal funds rate
% Sample period: 1970Q1 to 2011Q1
% Identification of 3 shocks, monetary policy shock, aggregate demand and 
% aggregate supply shock using sign restrictions.
% Algorithm follows Andrew Binning (2013)'s implementation of RWZ (2010)
% =======================================================================
% Willi Mutschler, December 2017
% willi@mutschler.eu
% =======================================================================

clearvars; clc;close all;

ENDO = xlsread('../data/data.xlsx','MonPolData'); % Load data
vnames   = {'GDP level','CPI','Federal Funds Rate'}; % Names for variables
epsnames = {'Aggregate Demand','Aggregate Supply','Monetary Policy'}; % Names for shocks
IRFcumsum = [1 1 0];      % Whether IRFs should be cumulated or not. [Kx1]
nlag = 4;                 % Number of lags. [Scalar]
const = 1;                % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend
EstMethodVAR = 'ML';      % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
StartValueMethod = 2;     % 0: Use identity matrix, 1: use square root, 2: use cholesky as starting value
nsteps = 21;              % Horizon for impulse responses
ndraws = 1000;            % Number of solutions that are supposed to match the sign restrictions
bupp = 0.84;              % Upper error band for impulse response (percentile)
blow = 0.16;              % Lower error band for impulse response (percentile)
it_print = 10;            % Print on the screen every "it_print"-th iteration

% Sign restrictions on impact
Rsign  = [+1  +1   -1;
          +1  -1   -1;
          +1  nan  +1;
         ];

% Estimate reduced-form model
if strcmp(EstMethodVAR,'OLS') % ordinary least squares
    VAR = VARReducedForm(ENDO,nlag,const);
    SIGMAUHAT = VAR.SIGu;
elseif strcmp(EstMethodVAR,'ML') % maximum likelihood
    VAR = VARReducedForm(ENDO,nlag,const);
    SIGMAUHAT = VAR.SIGMLu;
end
nvars = VAR.nvar;                           % Number of variables. [Scalar]
nobs = VAR.nobse;                           % Number of observations used in estimation.  [Scalar]
Acomp = VAR.Acomp;                          % Companion matrix of VAR representation. [nvars*nlags x nvars*nlags]

% Estimate structural model
I_K     = eye(nvars);
J       = [I_K zeros(nvars,nvars*(nlag-1))];

if StartValueMethod == 0
    B0inv_0 = eye(VAR.nvar); % Use identity matrix as starting value
elseif StartValueMethod == 1
    B0inv_0 = SIGMAUHAT^.5; % Use square root of vcov of reduced form as starting value
elseif StartValueMethod == 2
    B0inv_0 = chol(SIGMAUHAT,'lower'); % Use Cholesky decomposition of vcov of reduced form
end

% Initialize output array for accepted impulse resonse functions
% 1st dimension = variable
% 2nd dimension = shock
% 3rd dimension = horizon
% 4th dimension = draw
IRFvals = zeros(nvars,nvars,nsteps,ndraws); 

tic;
waitb = waitbar(0,'Number of accepted draws');
accepted_draws = 1;
discarded_draws = 0;
while accepted_draws < ndraws+1
    % Draw an orthogonal matrix Q from set set of all orthogonal matrices
    W = randn(nvars,nvars);
    [Q,R] = qr(W);
    % Normalization: make sure that diagonal of upper triangular R is positive
    for i = 1:nvars
        if R(i,i)<0           % if value on diagonal negative
            Q(:,i) = -Q(:,i); % reverse all signs in the i-th column of Q
        end
    end
    B0inv = B0inv_0*Q; % New draw for B0inv    
    
    % Check sign restrictions
    for jj = 1:nvars
        iota = I_K(:,jj);
        chk = B0inv*iota;                 % Get column of impact matrix
        THETA = zeros(nvars*nlag,nsteps); % initialize impulse response function used in companion form
        THETA(1:nvars,1) = B0inv*iota;    % initial impact of shock jj
        
        sr_index = ~isnan(Rsign(:,jj));   
        tmp = sign(chk(sr_index)) - Rsign(sr_index,jj); % Check restrictions
        % if sign restrictions are not fullfilled, return to beginning of while loop
        if any(tmp~=0)                
            jj = 0;
            discarded_draws = discarded_draws + 1; % Count discarded draws
            break 
        end
        % Else compute IRF for horizon up to nsteps given the companion form
        for ii = 2:nsteps
            THETA(:,ii) = Acomp*THETA(:,ii-1);                
        end
        % Save IRF of shock j into output structure, we only need first K rows of companion form
        IRFvals(:,jj,:,accepted_draws) = THETA(1:nvars,:);
    end
    
    if jj == nvars
        accepted_draws = accepted_draws + 1; % Count accepted draws
        if mod(accepted_draws,it_print) == 0 % Update waitbar every it_print
            waitbar(accepted_draws/ndraws);
        end
    end
end
close(waitb);
fprintf('Number of accepted draws: %d\n',accepted_draws-1);
fprintf('Number of discarded draws: %d\n',discarded_draws);
fprintf('Number of total draws: %d\n',accepted_draws-1+discarded_draws);
toc;

% Display ALL impulse response functions that satisfy restrictions
figure('Name','Impulse Response Functions');
countplots = 1;
for ishock = 1:nvars
    for ivar = 1:nvars
        subplot(nvars,nvars,countplots);
        IRFS = squeeze(IRFvals(ivar,ishock,1:nsteps,:));
        if IRFcumsum(ivar) == 1
            IRFS = cumsum(IRFS,1);
        end
        plot(0:1:nsteps-1, IRFS);
        title(vnames{ivar})
        ylabel([epsnames{ishock}, ' Shock'])
        countplots = countplots + 1;
    end
end

% Display median response of impulse response functions that satisfy restrictions
qstring = sprintf('Do you want to display median responses?\nNote that these:\n(1) lack structural economic interpretation\n      and\n(2) are not a useful statistical summary');
button = questdlg(qstring,'Median Responses','Yes, I don''t care','No','No');

if strcmp(button,'Yes, I don''t care')
    figure('Name','Median Response And Percentiles');
    countplots = 1;
    for ishock = 1:nvars
        for ivar = 1:nvars
            subplot(nvars,nvars,countplots);
            IRFS = squeeze(IRFvals(ivar,ishock,1:nsteps,:));
            if IRFcumsum(ivar) == 1
                IRFS = cumsum(IRFS,1);
            end        
            IRFSmed = median(IRFS,2);
            IRFSup  = prctile(IRFS,bupp*100,2);
            IRFSlow  = prctile(IRFS,blow*100,2);
            plot(0:1:nsteps-1, IRFSmed ,'b','linewidth',2);
            hold on;
            plot(0:1:nsteps-1, [IRFSlow IRFSup] ,'--r');
            grid;
            title(vnames{ivar})
            ylabel([epsnames{ishock}, ' Shock'])
            countplots = countplots + 1;
        end
    end
end