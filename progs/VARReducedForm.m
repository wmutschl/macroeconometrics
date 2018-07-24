function VAR = VARReducedForm(ENDO,nlag,const,dispestim)
% =======================================================================
% Perform vector autogressive (VAR) estimation with OLS
% =======================================================================
% VAR = VARReducedForm(ENDO,nlag,const)
% Representation: y_t = [c d A_1 ... A_p] [1 t y_{t-1}' ... y_{t-p}']' + u(t) = A Z_{t-1} + u_t
% -----------------------------------------------------------------------
% INPUT
%	- ENDO: an (nobs x nvar) matrix of endogenous variables,    
%           nobs: number of obserbations and nvar is number of variables
%	- nlag: lag length
%	- const: 0 no constant; 1 constant; 2 constant and linear trend
% -----------------------------------------------------------------------
% OUTPUT
%   - VAR: structure including VAR estimation results, see below for
%   details
% -----------------------------------------------------------------------
% CALLS
%   - OLSmodel.m: builtin function (see below) to estimate regression models with ols
%   - DisplayEstimation.m: builtin function (see below) to display estimation coefficients
% =======================================================================
% Willi Mutschler, November 2017
% willi@mutschler.eu
% Note: this code is a modified version of of the vare.m function of James
% P. LeSage and the VARmodel.m function of Cesa-Bianchi

if nargin < 4
    dispestim = 1;
end

%% Get some parameters
[nobs, nvar]  = size(ENDO);
nobse         = nobs - nlag;    % Effective sample size
ncoeff        = nvar*nlag;      % number of coefficients in A matrix
ntotcoeff     = ncoeff + const; % total number of coefficients to estimate

%% Create independent vector and lagged dependent matrix
ENDOt=transpose(ENDO);
% Y is matrixed of lagged endogenous variables
Y=ENDOt(:,nlag:nobs);
for i=1:nlag-1
 	Y=[Y; ENDOt(:,nlag-i:nobs-i)];
end
% add deterministic terms if any
if const == 0
    Z = Y(:,1:nobs-nlag);
elseif const == 1
    Z=[ones(1,nobs-nlag); Y(:,1:nobs-nlag)];
elseif const == 2
    Z=[ones(1,nobs-nlag); (nlag+1):nobs; Y(:,1:nobs-nlag)];
end
Y=Y(:,2:nobs-nlag+1);

%% OLS estimation equation by equation
for j=1:nvar
    y = Y(j,:)';
    x = Z';
    OLSout = OLSmodel(y,x);
    aux = sprintf('eq%d',j); % this creates strings 'eq1' 'eq2' 'eq3' ... 
                             %you can use this with (aux) to evaluate them 
                             % in order to name objects and stuff
    % compute t-probs
    tstat = OLSout.tstat;
    tout = tpdf(tstat,nobse-ncoeff);    
    VAR.(aux).beta  = OLSout.beta; % bhats
    VAR.(aux).tstat = OLSout.tstat; % t-stats
    VAR.(aux).tprob = tout;        % t-probs    
    VAR.(aux).resid = OLSout.resid;% resids 
    VAR.(aux).yhat  = OLSout.yhat; % yhats
    VAR.(aux).y     = y;           % actual y
    VAR.(aux).rsqr  = OLSout.rsqr; % r-squared
    VAR.(aux).rbar  = OLSout.rbar; % r-adjusted
    VAR.(aux).sige  = OLSout.sige; % standard error
end 

%% Compute the matrix of coefficients & Covariance Matrix
A = Y*Z'/(Z*Z');
U = Y-A*Z;
SIGMA = (1/(nobse-ntotcoeff))*(U*U'); % adjusted for # of estimated coeff per equation
% Companion matrix max eigenvalue
Acomp = [A(1:nvar,1+const:nvar*nlag+const); eye(nvar*(nlag-1)) zeros(nvar*(nlag-1),nvar)];
maxEig = max(abs(eig(Acomp)));
%% Save into structure
VAR.A = A(1:nvar,:);
VAR.Acomp = Acomp;
VAR.maxEig = maxEig;
VAR.SIGu = SIGMA(1:nvar,1:nvar);
VAR.SIGMLu = (nobse-ntotcoeff)/nobse*VAR.SIGu; % Maximum Likelihood COV Matrix
VAR.residuals = U;
VAR.Z = Z;
VAR.Y = Y;
VAR.ENDO = ENDO;
VAR.nlag = nlag;
VAR.const = const;
VAR.nobs  = nobs;
VAR.nobse = nobse;
VAR.nvar  = nvar;
VAR.nlag  = nlag;
VAR.ncoeff = ncoeff;
VAR.ntotcoeff = ntotcoeff;
VAR.const     = const;
if dispestim
    DisplayEstimation(VAR);
end

%% Additional Functions
function DisplayEstimation(V)
    Ahat = V.A;
    SigmahatU = V.SIGu;
    SigmatildeU = V.SIGMLu;
    cons = V.const;
    p = V.nlag;
    K = size(Ahat,1);
    % Display estimation results
    if cons == 0
        estimtable = table([]);
    elseif cons == 1
        nuhat = Ahat(:,1);
        estimtable = table(nuhat);
    elseif cons == 2
        nuhat = Ahat(:,1);
        timehat = Ahat(:,2);
        estimtable = table(nuhat,timehat);
    end
    ntrend = size(estimtable,2);
    Ai = reshape(Ahat(:,(1+ntrend):end),[K,K,p]);
    for ii = 1:p
        estimtable = [estimtable table(Ai(:,:,ii),'VariableNames', {sprintf('Ahat%d',ii)})];
    end
    disp(estimtable);
    disp([table(SigmahatU) table(SigmatildeU)]);
end % DisplayEstimation End

function OLS = OLSmodel(y,x,meth)
    % =======================================================================
    % OLS regression
    % =======================================================================
    % OLS = OLSmodel(y,x)
    % -----------------------------------------------------------------------
    % INPUT
    %	- y: dependent variable vector    (nobs x 1)
    %	- x: independent variables matrix (nobs x nvar)
    % -----------------------------------------------------------------------
    % OUPUT
    %	- OLS: structure including OLS estimation results
    % =======================================================================
    % Based on OLSmodel.m from Ambrogio Cesa Bianchi and olse.m from James P.
    % LeSage and fn_ols.m from Tao Tzha (Dynare implemenation).
    if nargin <3
        meth = 0; % Use SVD Decomposition
    end
    [T, K] = size(x); 

    if meth == 0 % Compute inv(X'X) with SVD
        [u d v]=svd(x,0);
        vd=v.*(ones(size(v,2),1)*diag(d)');
        dinv = 1./diag(d);    % inv(diag(d))
        vdinv=v.*(ones(size(v,2),1)*dinv');
        xtxinv = vdinv*vdinv';
        uy = u'*y;
        xty = vd*uy;
        beta = xtxinv*xty;
        yhat =  u*uy;    
    else
        if T < 10000 %Compute inv(X'X) with QR
            [~, r] = qr(x,0); 
            xtxinv = (r'*r)\eye(K);
        else %Compute inv(X'X) with Matlabs builtin function
            xtxinv = (x'*x)\eye(K);
        end
        beta = xtxinv*(x'*y);
        yhat = x*beta;    
    end
    resid = y - yhat;
    sigu = resid'*resid;
    sige = sigu/(T-K);
    tmp = (sige)*(diag(xtxinv));
    sigb=sqrt(tmp);
    tcrit=-tinv(0.025,T);
    bint=[beta-tcrit.*sigb, beta+tcrit.*sigb];
    tsta = beta./(sigb);

    ym = y - mean(y);
    rsqr1 = sigu;
    rsqr2 = ym'*ym;
    rsqr = 1.0 - rsqr1/rsqr2;
    rsqr1 = rsqr1/(T-K);
    rsqr2 = rsqr2/(T-1.0);
    if rsqr2 ~= 0
        rbar = 1 - (rsqr1/rsqr2);
    else
        rbar = rsqr;
    end
    ediff = resid(2:T) - resid(1:T-1);
    dw = (ediff'*ediff)/sigu; % durbin-watson


    OLS.beta = beta;
    OLS.yhat = yhat;
    OLS.resid = resid;
    OLS.sige = sige;
    OLS.bstd = sigb;
    OLS.bint=bint;
    OLS.tstat = tsta;
    OLS.rsqr = rsqr;
    OLS.rbar = rbar;
    OLS.dw = dw;
    OLS.y = y;
    OLS.x = x;
    OLS.nobs = T;
    OLS.nvar = K;

end%OLSmodel end


end %Main Function end