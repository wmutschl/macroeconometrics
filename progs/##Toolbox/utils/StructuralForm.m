function [StructMod,opt] = StructuralForm(ReducedForm,opt)
% =======================================================================
% Estimates the structural model (SVAR or SVECM) provided an estimate of 
% the reduced-form model. Algorithms for deriving the structural model:
%    - Cholesky decomposition for recursive models
%    - Numerical optimization using Matlab's fsolve
%    - QR algorithm of Rubio-Ramirez,Waggoner,Zha (2010) or Binning (2013)
%    - Scoring algorithm of Amisano and Giannini (1997)
% The normalization rule imposes taht the sign of the diagonal elements
% of B0inv are positive.
% If the model is underidentified using only exclusion restrictions, sign 
% restrictions are used to narrow down the set of valid models.
% =======================================================================
% [StructMod,opt] = StructuralForm(ReducedForm,opt)
% -----------------------------------------------------------------------
% INPUTS
%	- ReducedForm : structure including reduced form estimation results, see load reduced form section below
%   - opt  : structure with options, see load load options section below
% -----------------------------------------------------------------------
% OUTPUTS
%   - StructMod : structure including structural estimation with all fields of ReducedForm and the following fields:
%       - B0inv_point : Estimate of short-run impact matrix. [nvars x nvars]
%       - LRMat_point : Estimate of long-run impact matrix. [nvars x nvars]
% -----------------------------------------------------------------------
% CALLS
%   - f_SMOD
%   - findQs
%   - generateDraw
%   - commutation
%   - vec
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

% Load options
model = opt.model;                       % Reduced form model, values are 1 for VAR or 2 for VECM. [scalar]
algorithm = opt.algorithm;               % Algorithm used for estimation, values are 'Cholesky', 'Optimization', 'RWZ', 'Scoring'. [String] 
StartValueMethod = opt.StartValueMethod; % Starting values of algorithms, values are 1 or 2. [Scalar]
                                         % - RWZ & Optimization: 1 square root or 2 Cholesky decomposition of covariance matrix of reduced-form residuals
                                         % - Scoring: 1 draw randomly, 2 fixed values (see StructuralForm.m)
nlag = opt.nlag;                         % Number of lags. [Scalar]

Rshort = opt.Rshort;                   % Matrix of short-run restrictions on impact matrix B_0^{-1}. [nvars x nvars]
Rshorth = opt.Rshorth;                   % Matrix of short-run restrictions on impact matrix at a finite horizon h, i.e. on (J*Acomp^h*J')*B_0^{-1}. [nvars x nvars]
Rlong = opt.Rlong;                       % Matrix of long-run restrictions. [nvars x nvars]
                                         % - for VAR on A(1)*B0inv = inv(eye(nvars)-A1hat-A2hat-...-Aphat))B0inv
                                         % - for VECM on Upsylon*B0inv = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))*B0inv,  
Rsign = opt.Rsign;                       % Matrix of sign restrictions. [nvars x nvars]
Rsignh = opt.Rsignh;                       % Matrix of sign restrictions. [nvars x nvars]

% Load reduced-form estimation
nvars = size(ReducedForm.ENDO,2);        % Number of variables. [Scalar]
nobs = ReducedForm.nobs;                 % Number of observations used in estimation.  [Scalar]
SIGMAUHAT = ReducedForm.SIGMAUHAT;       % Estimated covariance matrix of reduced-form residuals. [nvars x nvars]
Acomp = ReducedForm.Acomp;               % Companion matrix of VAR representation. [nvars*nlags x nvars*nlags]
if model == 2 % for VECM model
    alph = ReducedForm.paramVals.A;      % Loading matrix. [nvars x coint_r]
    betta = ReducedForm.paramVals.B;     % Matrix of cointegration vectors. [nvars x coint_r]
    paramVals = ReducedForm.paramVals;   % Structure of parameter estimates of VECM model
end

% Precomputations
selSR  = find(isnan(Rshort)==0); % Index for short-run restrictions on impact
selLR  = find(isnan(Rlong)==0);  % Index for long-run restrictions
selSRh = {}; idxSRh = [];
for i = 1:size(Rshorth,2)
    if isempty(Rshorth{i}) == 0
        idxSRh = [idxSRh i];
        selSRh{i} = find(isnan(Rshorth{i})==0); % Index for short-run restrictions on impact
    end
end
selSign = find(isnan(Rsign)==0); % Index for sign restrictions on impact
selSignh = {}; idxSignh = [];
for i = 1:size(Rsignh,2)
    if isempty(Rsignh{i}) == 0
        idxSignh = [idxSignh i];
        selSignh{i} = find(isnan(Rsignh{i})==0); % Index for sign restrictions
    end
end

I_K    = eye(nvars);
I_KK   = eye(nvars^2);
J = [I_K zeros(nvars,nvars*(nlag-1))];
    
if model == 1 %VAR model
    A1inv_big = inv(eye(size(Acomp,1))-Acomp); % from the companion
    LRMat = A1inv_big(1:nvars,1:nvars);        % total impact matrix A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)
elseif model == 2 %VECM model
    GAM = I_K;
    if nlag >0
        for i = 1:(nlag-1)
            GAM = GAM - paramVals.(sprintf('B%d',i)); % compute mean lag matrix
        end
    end
    betta_o = null(betta');	%orthogonal complement of beta
    alph_o = null(alph');	%orthogonal complemt of alpha
    LRMat = betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)); %total impact matrix Upsylon = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))
end

StructMod = ReducedForm;

%% Exactly or over-identified models without sign restrictions
if (isempty(selSign) == 1) && (isempty(selSignh) == 1)
    StructMod.SignRestr = 0;
    switch algorithm

        case 'Cholesky'  % Cholesky decomposition for recursive structural models
            if (isempty(selSR) == 0) && (isempty(selSRh) == 1) && (isempty(selLR) == 1)
                B0inv = chol(SIGMAUHAT,'lower'); % Only recursive short-run restrictions on impact
            elseif (isempty(selSR) == 1) && (isempty(selSRh) == 1) && (isempty(selLR) == 0)
                B0inv = LRMat\chol(LRMat*SIGMAUHAT*LRMat','lower'); % Only recursive long-run restrictions  
            else
                error('In order to use Cholesky, please specify recursive structure for only short-run at impact or only long-run, not for both. Or use a different algorithm for nonrecursive structures.');
            end

        case 'Optimization'  % Numerical Optimization using fsolve
            % Options for fsolve
            TolX = 1e-4;            % Termination tolerance on the current point
            TolFun = 1e-9;          % Termination tolerance on the function value
            MaxFunEvals = 1e+50000; % Maximum number of function evaluations allowed
            MaxIter = 1000;         % Maximum numberof iterations allowed
            NumIter = 3;            % Number of iterations on starting values
            OptimAlgorithm = 'trust-region-dogleg'; % Algorithm used in fsolve
            options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
            if StartValueMethod == 1
                B0inv_0= SIGMAUHAT^.5; % Use square root of vcov of reduced form as starting value
            elseif StartValueMethod == 2
                B0inv_0 = chol(SIGMAUHAT,'lower'); % Use Cholesky decomposition of vcov of reduced form
            end        
            % Call optimization routine fsolve to minimize f_SMOD.m
            [B0inv,fval,exitflag,output]=fsolve('f_SMOD',B0inv_0,options,SIGMAUHAT,LRMat,J,Acomp,Rshort,Rshorth,Rlong,selSR,selSRh,selLR,idxSRh);
            if strcmp(output.algorithm,'levenberg-marquardt') % Catch possible change of algorithms to surpress warnings
                options.Algorithm = output.algorithm;            
            end
            for r = 1:NumIter % Iterate on starting values
                [B0inv,fval,exitflag,output]=fsolve('f_SMOD',B0inv,options,SIGMAUHAT,LRMat,J,Acomp,Rshort,Rshorth,Rlong,selSR,selSRh,selLR,idxSRh);
            end

        case 'RWZ' % QR algorithm of Rubio-Ramirez,Waggoner and Zha (2010), for implementation details see Binning (2013)        
            if sum(Rshort(selSR)) ~= 0 || sum(Rlong(selLR)) ~= 0 || sum(sum(cell2mat(Rshorth),'omitnan')) ~= 0
                error('RWZ algorithm does not support nonzero restrictions')
            end
            if StartValueMethod == 1
                B0inv_0 = chol(SIGMAUHAT,'lower'); % Use Cholesky decomposition of vcov of reduced form
            elseif StartValueMethod == 2
                B0inv_0= SIGMAUHAT^.5; % Use square root of vcov of reduced form as starting value
            end        
            % Set up linear restrictions        
            f = [Rshort;cell2mat(transpose(Rshorth));Rlong]; % Translate transformation function into a matrix of zeros and ones
            f(isnan(f)) = 1; % set NaN to 1        
            [Q,index,flag] = findQs(nvars,f); % Find Q matrices that correspond to the linear restritions on the columns of the transformed space
                                              % Note that the order of variables is changed according to index

            B0inv = generateDraw(B0inv_0,nvars); % Generate draw that is consistent with the shock covariance matrix
            L0 = B0inv; % Initial short-run impact matrix
            Lh = [];
            for h = idxSRh
                Lh = [Lh;(J*Acomp^h*J')*L0];
            end
            Linf = LRMat*B0inv; % Initial long-run impact matrix
            F = [L0;Lh;Linf]; % Initial impact matrix
            P = zeros(nvars,nvars); % Initialize rotation matrix
            % Find orthogonal rotation matrix P, such that F*P = f  and P*P'=I
            for ii = 1:nvars
                if ii == 1
                    Qtilde = Q{ii}*F;
                else
                    Qtilde = [Q{ii}*F;P'];
                end
                [QQ,RR] = qr(Qtilde');
                P_temp = QQ(:,end);
                P(:,ii) = P_temp;    
            end
            B0inv = L0*P(:,index); % Variables are reordered to original ordering

        case 'Scoring' % Scoring algorithm as in Amisano and Gianini (1997), Vlaar (2004) or Breitung, Brueggemann and Luetkepohl (2004)
            % Options for scoring algorithm
            eps1_tol = 1E-6;	   % Tolerance of relative change in parameters 
            eps2_tol = 1E-10; 	   % Tolerance for relative change in log likelihood
            maxIterations = 500;   % Maximum number of iterations
            fixStart = .1;         % fixed starting value if StartValueMethod = 2         

            % Compute restriction matrices in implicit form R*vec(B0inv) = r
            R_SR  = I_KK(selSR,:);
            r_SR  = Rshort(selSR);
            R_SRh = [];
            r_SRh = [];
            for h = idxSRh            
                rshorth = Rshorth{h};
                R_SRh = [R_SRh;I_KK(selSRh{h},:)*kron(I_K,J*Acomp^h*J')];
                r_SRh = [r_SRh; rshorth(selSRh{h})];
            end
            R_LR  = I_KK(selLR,:);
            r_LR  = Rlong(selLR);
            if isempty(R_LR)
                R_LR = [];
                r_LR = [];
            else	
                R_LR = R_LR*kron(I_K,LRMat); % Restate long-run restrictions only on B0inv
            end
            if isempty(R_SR) == 0
                R = [R_LR;R_SR];
                r = [r_LR;r_SR];
            else
                R = R_LR;
                r = r_LR;
            end
            R = [R;R_SRh];
            r = [r;r_SRh];
            % Compute restriction matrices in explicit form vec(B0inv) = null(R)*gam + R'*inv(R*R')*r = R_o*gam + r_o where gam denotes the free parameters
            R_o = null(R); % Compute null space of R, i.e. R*R_o = 0
            r_o = pinv(R)*r;

            % Algorithm for restrictions for AB-model: [vec(A);vec(B)] = [R_A 0; 0 R_B]*[gam_A;gam_B] + [r_A r_B],
            % see e.g. Breitung, Brueggemann and Luetkepohl (2004). Note that we have a B-model, therefore we restrict 
            % the A matrix to the identity matrix. Note further that in our setting B0inv = inv(A)*B.

            % Set up [vec(A);vec(B)] = [R_A 0; 0 R_B]*[gam_A;gam_B] + [r_A r_B]
            R_A  = zeros(nvars^2,size(R_o,2));
            R_B  = R_o;
            R_AB = [R_A;R_B];
            r_A  = I_K(:);
            r_B  = r_o;
            r_AB = [r_A;r_B];

            % Count restrictions for identification
            nfreeall = size(R_o,2);
            if isempty(R_LR) == 0
                nResLRMat  = rank(R_LR);
            else
                nResLRMat = 0;
            end
            if isempty(R_SR) == 0
                nResB0inv = rank(R_SR);
            else
                nResB0inv = 0;
            end
            if isempty(R_SRh) == 0
                nResSRh = rank(R_SRh);
            else
                nResSRh = 0;
            end

            % Starting Values for algorithm
            if StartValueMethod == 1
                gam_AB = randn(size(R_AB,2),1)*.1; %random draw
            elseif StartValueMethod == 2
                gam_AB = ones(size(R_AB,2),1)*fixStart; %fixed values
            end

            % Check identification
            ident = 0; % Set identification flag to 0
            if (nResB0inv+nResLRMat+nResSRh) < (nvars*(nvars-1)/2)
                DegreeofOverIdent = 0;
            else 	 
                vecAB = R_AB*gam_AB+r_AB;
                At = reshape(vecAB(1:nvars^2),nvars,nvars)';
                Bt = reshape(vecAB(nvars^2+1:2*nvars^2),nvars,nvars)';
                Com = commutation(nvars,nvars);
                Btinv = Bt\I_K;
                V   = [(I_KK+Com)*kron(((At\Bt))',Btinv)  -1*(I_KK+Com)*kron(I_K,Btinv)];
                DegreeofOverIdent = nvars*(nvars+1)/2 - nfreeall;

                if isempty(V)==0 && isempty(R_AB)==0
                    if sum(eig((V*R_AB)'*(V*R_AB)) < 1E-10) == 0
                        ident = 1;
                    end
                end	  	
            end

            if ident == 1
                % Initialization
                maxls   = 1; % maximum for likelihood step
                vecAB   = R_AB*gam_AB+r_AB; % Restrictions in explicit form
                A   	= reshape(vecAB(1:nvars^2,1),nvars,nvars); % current A matrix
                B       = reshape(vecAB(nvars^2+1:2*nvars^2,1),nvars,nvars); % current B matrix
                Com 	= commutation(nvars,nvars); % commutation matrix
                noConvergence = 0; % set convergence flag to everything is fine
                gam_AB_old   = gam_AB; % initialize gamma vector
                eps1 	= 100; % Set tolerance of relative change in parameters high
                eps2 	= 100; % Set Tolerance for relative change in log likelihood high
                AinvB=A\B; 
                SigmaUt = AinvB*AinvB';
                BinvA = B\A;
                llold = nobs/2*log(det(BinvA)^2)-nobs/2*sum(diag((BinvA'*BinvA*SigmaUt))); % Initialize likelihood

                % Scoring algorithm
                i = 0;
                while (isempty(find(eps1 > eps1_tol))==0 || (eps2 > eps2_tol)) && i<maxIterations
                    % Compute information matrix and score vectors
                    BinvA = B\A;
                    Btinv = transpose(B)\eye(size(B,2));
                    InfoMat_vecAB = nobs*[ kron(inv(BinvA),Btinv) ; -1*kron(I_K,Btinv) ]*(I_KK+ Com)*[ kron(inv(BinvA'),inv(B)) -1*kron(I_K,inv(B))]; % Information matrix for [vec(A);vec(B)]
                    InfoMat_gam_AB = R_AB'*InfoMat_vecAB*R_AB; % Information matrix for [gamA;gamB]
                    invBinvAt = inv(BinvA)';
                    Score_vecBinvA  = nobs*(invBinvAt(:)'-transpose(BinvA(:))*kron(SIGMAUHAT,I_K));
                    Score_vecAB     = Score_vecBinvA * [kron(I_K,inv(B))  -1*kron(A'*Btinv,inv(B))];
                    Score_gam_AB    = Score_vecAB*R_AB;  
                    InfoABinv_times_Score_gamAB = InfoMat_gam_AB\transpose(Score_gam_AB);                
                    length = max(abs(InfoABinv_times_Score_gamAB));
                    if length > maxls % Check whether step size is largen than one, if so then flip                
                        lambda = maxls/length;
                    else
                        lambda = 1;
                    end
                    % Iterate on free parameters and compute A, B matrices as well as new likelihood value
                    gam_AB = gam_AB_old +lambda*InfoABinv_times_Score_gamAB; % Iterate
                    vecAB = R_AB*gam_AB+r_AB;
                    A = reshape(vecAB(1:nvars^2),nvars,nvars);
                    B = reshape(vecAB(nvars^2+1:2*nvars^2),nvars,nvars);   
                    BinvA = B\A;
                    AinvB=A\B;
                    SigmaUt = AinvB*AinvB';
                    llnew  = nobs/2*log(det(BinvA)^2)-nobs/2*sum(diag((BinvA'*BinvA*SigmaUt)));

                    eps2 = abs((llnew-llold)/llold); % Compute relative change in likelihood
                    eps1 = abs((gam_AB-gam_AB_old)./gam_AB_old); % Compute relative change in free parameters
                    gam_AB_old = gam_AB; 
                    llold  = llnew;
                    i = i +1;  
                end
                B0inv = A\B; % This is due to our setting

                if isempty(find(eps1>eps1_tol))==0 || (eps2>eps2_tol);
                  noConvergence = 1; % Something went wrong
                end
            end % if ident==1 end

    end % switch algorithm end

    % Normalize sign of B0inv such that diagonal elements are positive
    if sum(diag(B0inv)<0) ~= 0
        x = diag(B0inv)<0;
        B0inv(:,find(x==1)) = -1*B0inv(:,find(x==1));
    end
    StructMod.B0inv_point = B0inv;
    StructMod.LRMat_point = LRMat*B0inv;

else
%% Models with sign restrictions use QR algorithm of Binning (2013)
    StructMod.SignRestr = 1;
    if sum(Rshort(selSR)) ~= 0 || sum(Rlong(selLR)) ~= 0 || sum(sum(cell2mat(Rshorth),'omitnan')) ~= 0
        error('RWZ algorithm does not support nonzero restrictions')
    end
    if StartValueMethod == 1
        B0inv_0 = chol(SIGMAUHAT,'lower'); % Use Cholesky decomposition of vcov of reduced form
    elseif StartValueMethod == 2
        B0inv_0= SIGMAUHAT^.5; % Use square root of vcov of reduced form as starting value
    end        
    % Set up linear restrictions        
    f = [Rshort;cell2mat(transpose(Rshorth));Rlong]; % Translate transformation function into a matrix of zeros and ones
    f(isnan(f)) = 1; % set NaN to 1        
    [Q,index,flag] = findQs(nvars,f); % Find Q matrices that correspond to the linear restritions on the columns of the transformed space
                                      % Note that the order of variables is changed according to index
    counter = 1;         
    reverseStr = '';
    draws = opt.ndraws; % Number of draws
    T = opt.nsteps; % length of impulse response function.
    R = zeros(nvars,T,draws,nvars); % Contains impulse resonse functions
    % 1st dimension = variable
    % 2nd dimension = time
    % 3rd dimension = draw
    % 4th dimension = shock
    while counter < draws+1
        B0inv = generateDraw(B0inv_0,nvars); % Generate draw that is consistent with the shock covariance matrix
    
        L0 = B0inv; % Initial short-run impact matrix
        Lh = [];
        for h = idxSRh
            Lh = [Lh;(J*Acomp^h*J')*L0];
        end
        Linf = LRMat*B0inv; % Initial long-run impact matrix
        F = [L0;Lh;Linf]; % Initial impact matrix
        P = zeros(nvars,nvars); % Initialize rotation matrix
        % Find orthogonal rotation matrix P, such that F*P = f  and P*P'=I
        for ii = 1:nvars
            if ii == 1
                Qtilde = Q{ii}*F;
            else
                Qtilde = [Q{ii}*F;P'];
            end
            [QQ,RR] = qr(Qtilde');
            P_temp = QQ(:,end);
            P(:,ii) = P_temp;    
        end
        P = P(:,index);
        W = B0inv*P;

        for jj = 1:nvars
% Efficiency! skip for now...
%             if W(var_pos(jj),jj) < 0
%                 shock = -I_K(:,jj);
%             else
%                 shock = I_K(:,jj);
%             end
            shock = I_K(:,jj);
            V = zeros(nvars*nlag,T);

            V(1:nvars,1) = W*shock;

            chk = W*shock;
            sr_index = ~isnan(Rsign(:,jj));
            tmp = sign(chk(sr_index)) - Rsign(sr_index,jj);

            if any(tmp~=0)                
                jj = 0;
                break
            end

            for ii = 2:T
                V(:,ii) = Acomp*V(:,ii-1);                
            end
            R(:,:,counter,jj) = V(1:nvars,:);

        end
        if jj == nvars
            if counter < draws
                msg = sprintf('%d out of %d draws found that satisfy the sign and exclusion restrictions',counter,draws); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, size(msg,2));
            elseif counter == draws
                msg = sprintf('%d out of %d draws found that satisfy the sign and exclusion restrictions\n',counter,draws); fprintf([reverseStr, msg]); reverseStr = repmat(sprintf('\b'), 1, size(msg,2));
            end
            counter = counter + 1;            
        end

    end
    StructMod.IRFs = R;
end