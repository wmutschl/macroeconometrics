function VAR = ReducedForm(ENDO,nlags,VARopt)
% =======================================================================
% Perform vector autogressive (VAR) estimation with OLS 
% =======================================================================
% [VAR, VARopt] = VARmodel(ENDO,nlag,const,EXOG,nlag_ex)
% -----------------------------------------------------------------------
% INPUT
%	- ENDO: an (nobs x nvar) matrix of y-vectors
%	- nlag: lag length
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%	- const: 0 no constant; 1 constant; 2 constant and trend; 3 constant, 
%       trend, and trend^2 [dflt = 0]
%	- EXOG: optional matrix of variables (nobs x nvar_ex)
%	- nlag_ex: number of lags for exogeonus variables [dflt = 0]
% -----------------------------------------------------------------------
% OUTPUT
%   - VAR: structure including VAR estimation results
%   - VARopt: structure including VAR options (see VARoption)
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% Note: this code is a modified version of of the vare.m function of James 
% P. LeSage

% Representation -->  Y = Y(-1)*F' + u

% Note: compared to Eviews, there is a difference in the estimation of the 
% constant when lag is > 2. This is because Eviews initialize the trend
% with the number of lags (i.e., when lag=2, the trend is [2 3 ...T]), 
% while VARmakexy.m initialize the trend always with 1.

% I thank Jan Capek for spotting and addressing a compatibility issue with
% Matlab R2014a

% Adapted to fit notation of Kilian and Lütkepohl


%===============================================
DetermTrend = VARopt.DetermTrend;
[nvar,nobs] = size(ENDO);
nobse         = nobs - nlags; %actual observations used in estimation
ncoeff        = nvar*nlags;   %actual number of coefficients used in estimation
ntotcoeff     = ncoeff + DetermTrend; % total number of coefficients including deterministics trend variables
% Create independent vector and lagged dependent matrix
[Y, Z] = MakeYZ(ENDO,nlags,DetermTrend);

switch VARopt.EstimationMethod
    case 1 % OLS estimation
    %if nobs < 10000 % use svd or qr  
    %else
    Z_Zp_inv = Z*transpose(Z)\eye(ntotcoeff);
    %end;
    B = Y*Z'*Z_Zp_inv;
    Yhat = B*Z;
    Uhat = Y - Yhat;
    SigmaU = Uhat*transpose(Uhat)/(nobse-ntotcoeff);    
    SigmaB = diag(kron(Z_Zp_inv,SigmaU));
    tcrit=-tdis_inv(.025,nobse-ntotcoeff);
    Bintlow  = reshape(B(:)-tcrit.*sqrt(SigmaB),size(B));
    Binthigh = reshape(B(:)+tcrit.*sqrt(SigmaB),size(B));
    tstat = reshape(B(:)./sqrt(SigmaB),size(B));
    tout = tdis_prb(tstat,nobse-ntotcoeff);
end
% Companion Form
Acomp = [B(:,1+DetermTrend:nvar*nlags+DetermTrend); eye(nvar*(nlags-1)) zeros(nvar*(nlags-1),nvar)];
maxEig = max(abs(eig(Acomp)));

% Create Structure
VAR.ENDO = ENDO;
VAR.DetermTrend = DetermTrend;
VAR.nobse       = nobse;
VAR.nvar        = nvar;
VAR.nlags       = nlags;
VAR.ncoeff      = ncoeff;
VAR.ntotcoeff   = ntotcoeff;
VAR.Y           = Y;
VAR.Z           = Z;
VAR.B           = B;
VAR.SigmaU      = SigmaU;
VAR.Uhat        = Uhat;
VAR.Z           = Z;
VAR.Y           = Y;
VAR.Acomp       = Acomp;
VAR.maxEig      = maxEig;
% Some additional equation by equation results
for j=1:nvar;        
    VAR.(sprintf('eq%d',j)).beta  = B(j,:);      % bhats    
    VAR.(sprintf('eq%d',j)).tstat = tstat(j,:);  % t-stats
    VAR.(sprintf('eq%d',j)).tprob = tout(j,:);   % t-probs
    VAR.(sprintf('eq%d',j)).Bintlow = Bintlow(j,:);   % CI low
    VAR.(sprintf('eq%d',j)).Binthigh = Binthigh(j,:);   % CI high
    VAR.(sprintf('eq%d',j)).u     = Uhat(j,:);   % resids 
    VAR.(sprintf('eq%d',j)).yhat  = Yhat(j,:);   % yhats
    VAR.(sprintf('eq%d',j)).y     = Y(j,:);      % actual y
    VAR.(sprintf('eq%d',j)).sigu  = SigmaU(j,j); % standard error
end 

if VAR.maxEig>=1
    error('CheckStabilityCompForm:Stability','The VAR is not stationary!')
end;

%% Initialize other results
%===============================================
VAR.invA = [];  % inverse of teh A matrix (need identification: see VARir/VARfevd)
VAR.S    = [];  % Orthonormal matrix (need identification: see SR)
