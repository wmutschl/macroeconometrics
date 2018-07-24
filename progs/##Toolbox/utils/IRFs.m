function IRFpoint  = IRFs(StructMod,opt)
% =======================================================================
% Compute IRFs for a VAR model estimated with VARmodel. Three
% identification schemes can be specified: zero short-run restrictions,
% zero long run restrictions, and sign restrictions
% =======================================================================
% [IRF, VAR] = VARir(VAR,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - VAR: structure, result of VARmodel function
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% ----------------------------------------------------------------------- 
% OUTPUT
%   - IRF(t,j,k): matrix with 't' steps, containing the IRF of 'j' variable 
%       to 'k' shock
%   - VAR: structure including VAR estimation results
%       * VAR.invA: : identified contemporaneous A matrix
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% Note. This code follows the notation as in the lecture notes available at
% https://sites.google.com/site/ambropo/MatlabCodes

% Specifically if Y, u, and e are [NxT] matrices, the following is true:
% Reduced form VAR    -->  Y = F*Y(-1) + u    (F => Ft' from VARmodel function) 
% Structural form VAR --> AY = B*Y(-1) + e     

% Where:
%     F = invA*B  -->  B = A*F  -->  B = IRF.invA\transpose(VAR.Ft)
%     u = invA*e  -->  e = A*u  -->  e = IRF.invA\transpose(VAR.residuals);

% Impulse responses:
% IRF(1) = invA*e          where "e" is impulse in the code
% IRF(j) = H(j)*IRF(1)     where H(j)=H^j and for j=2,...,nsteps



%% Check inputs
%===============================================
%% Retrieve and initialize variables 
%===============================================
nsteps = opt.nsteps;
impact = opt.impact;
nvars  = size(StructMod.ENDO,2);
nlag   = opt.nlag;
B0inv  = StructMod.B0inv_point;
Acomp  = StructMod.Acomp;

IRFpoint = nan(nsteps,nvars,nvars);


%% Compute the impulse response
%===============================================
for mm=1:nvars

    % Initialize the impulse response vector
    response = zeros(nvars, nsteps);
    
    % Create the impulse vector
    impulse = zeros(nvars,1);
    
    % Set the size of the shock
    if impact==0
        impulse(mm,1) = 1; % one stdev shock
    elseif impact==1
        impulse(mm,1) = 1/B0inv(mm,mm); % unitary shock
    else
        error('Impact must be either 0 or 1');
    end

    % First period impulse response (=impulse vector)
    response(:,1) = B0inv*impulse;

    % Make it comparable with companion
    impulse_big  = [response(:,1)' zeros(1, nvars*(nlag-1))]';
    
    % Recursive computation of impulse response
    Acomp_eye = eye(size(Acomp,1)); 
    for kk = 2:nsteps
        Acomp_eye = Acomp * Acomp_eye; % this is the multiplier Fcomp^n
        response_big   = Acomp_eye * impulse_big;
        response(:,kk) = response_big(1:nvars);
    end
    IRFpoint(:,:,mm) = response';
end