function FEVDpoint = FEVDs(StructMod,opt)
% =======================================================================
% Compute FEVDs for a VAR model estimated with VARmodel. Three
% identification schemes can be specified: zero short-run restrictions,
% zero long run restrictions, and sign restrictions
% =======================================================================
% [FEVD, VAR] = VARfevd(VAR,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% OUTPUT
%	- FEVD(t,j,k): matrix with 't' steps, the FEVD due to 'j' shock for 
%       'k' variable
%   - VAR: structure including VAR estimation results
%       * VAR.invA: identified contemporaneous A matrix
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% I thank Dora Xia for pointing out a typo in the above description.


%% Check inputs
%===============================================
%% Retrieve and initialize variables 
%===============================================
nsteps = opt.nsteps;
nvars  = size(StructMod.ENDO,2);
SIGMAUHAT = StructMod.SIGMAUHAT;

MSE   = zeros(nvars,nvars,nsteps);
MSE_j = zeros(nvars,nvars,nsteps);
PSI   = zeros(nvars,nvars,nsteps);
FEVDpoint  = zeros(nsteps,nvars,nvars);
SE    = zeros(nsteps,nvars);

B0inv_point = StructMod.B0inv_point;
%% Compute the multipliers
%===============================================
Modjunk = StructMod;
Modjunk.SIGMAUHAT = eye(nvars);
[StructModjunk,opt] = StructuralForm(Modjunk,opt);
IRFjunk = IRFs(StructModjunk,opt);

% this loop is to save the multipliers for each period
for mm = 1:nvars
    PSI(:,mm,:) = reshape(IRFjunk(:,:,mm)',1,nvars,nsteps);
end


%% Calculate the contribution to the MSE for each shock (i.e, FEVD)
%===============================================

for mm = 1:nvars % loop for the shocks
    
    % Calculate Total Mean Squared Error
    MSE(:,:,1) = SIGMAUHAT;

    for kk = 2:nsteps;
       MSE(:,:,kk) = MSE(:,:,kk-1) + PSI(:,:,kk)*SIGMAUHAT*PSI(:,:,kk)';
    end;
    
    % Get the column of invA corresponding to the mm_th shock
    column = B0inv_point(:,mm);
    
    % Compute the mean square error
    MSE_j(:,:,1) = column*column';
    for kk = 2:nsteps
        MSE_j(:,:,kk) = MSE_j(:,:,kk-1) + PSI(:,:,kk)*(column*column')*PSI(:,:,kk)';   
    end

    % Compute the Forecast Error Covariance Decomposition
    FECD = MSE_j./MSE;
    
    % Select only the variance terms
    for nn = 1:nsteps
        for ii = 1:nvars
            FEVDpoint(nn,mm,ii) = FECD(ii,ii,nn);
            SE(nn,:) = sqrt( diag(MSE(:,:,nn))' );
        end
    end
end
