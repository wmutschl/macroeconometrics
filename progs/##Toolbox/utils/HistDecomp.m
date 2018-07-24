function HD = HistDecomp(StructMod,opt)
% =======================================================================
% Compute the historical decomposition of the time series in a VAR
% estimated with VARmodel and identified with VARir/VARfevd
% =======================================================================
% HD = VARhd(VAR)
% -----------------------------------------------------------------------
% INPUTS 
%   - VAR: structure, result of VARmodel -> VARir/VARfevd function
% -----------------------------------------------------------------------
% OUTPUT
%   - HD: structure including the historical decomposition
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com

% I thank Andrey Zubarev for finding a bug in the contribution of the 
% exogenous variables when nvar_ex~=0 and nlag_ex>0. 


%% Check inputs
%===============================================
B0inv = StructMod.B0inv_point;

%% Retrieve and initialize variables 
%===============================================
const    = opt.const;                     % constant and/or trends
nvars    = size(StructMod.ENDO,2);                     % number of endogenous variables
nvars_ex = size(StructMod.EXOG,2);                  % number of exogenous (excluding constant and trend)
nlag     = opt.nlag;                      % number of lags 
nlag_ex  = opt.nlag_ex;                   % number of lags of the exogenous 
nvarXeq  = nvars * nlag;                  % number of lagged endogenous per equation

nobs    = StructMod.nobs;                            % number of observations
Acomp   = StructMod.Acomp;                % Companion matrix
AMAT    = StructMod.AMATt';               % make comparable to notes
eps     = B0inv\transpose(StructMod.residuals);  % structural errors 
Y       = StructMod.Y;                          % left-hand side
X       = StructMod.X(:,1+const:nvarXeq+const); % right-hand side (no exogenous)



%% Compute historical decompositions
%===============================================

% Contribution of each shock
    invA_big = zeros(nvarXeq,nvars);
    invA_big(1:nvars,:) = B0inv;
    Icomp = [eye(nvars) zeros(nvars,(nlag-1)*nvars)];
    HDshock_big = zeros(nlag*nvars,nobs+1,nvars);
    HDshock = zeros(nvars,nobs+1,nvars);
    for j=1:nvars; % for each variable
        eps_big = zeros(nvars,nobs+1); % matrix of shocks conformable with companion
        eps_big(j,2:end) = eps(j,:);
        for i = 2:nobs+1
            HDshock_big(:,i,j) = invA_big*eps_big(:,i) + Acomp*HDshock_big(:,i-1,j);
            HDshock(:,i,j) =  Icomp*HDshock_big(:,i,j);
        end
    end
    
% Initial value
    HDinit_big = zeros(nlag*nvars,nobs+1);
    HDinit = zeros(nvars, nobs+1);
    HDinit_big(:,1) = X(1,:)';
    HDinit(:,1) = Icomp*HDinit_big(:,1);
    for i = 2:nobs+1
        HDinit_big(:,i) = Acomp*HDinit_big(:,i-1);
        HDinit(:,i) = Icomp *HDinit_big(:,i);
    end
    
% Constant
    HDconst_big = zeros(nlag*nvars,nobs+1);
    HDconst = zeros(nvars, nobs+1);
    CC = zeros(nlag*nvars,1);
    if const>0
        CC(1:nvars,:) = AMAT(:,1);
        for i = 2:nobs+1
            HDconst_big(:,i) = CC + Acomp*HDconst_big(:,i-1);
            HDconst(:,i) = Icomp * HDconst_big(:,i);
        end
    end
    
% Linear trend
    HDtrend_big = zeros(nlag*nvars,nobs+1);
    HDtrend = zeros(nvars, nobs+1);
    TT = zeros(nlag*nvars,1);
    if const>1;
        TT(1:nvars,:) = AMAT(:,2);
        for i = 2:nobs+1
            HDtrend_big(:,i) = TT*(i-1) + Acomp*HDtrend_big(:,i-1);
            HDtrend(:,i) = Icomp * HDtrend_big(:,i);
        end
    end
    
% Quadratic trend
    HDtrend2_big = zeros(nlag*nvars, nobs+1);
    HDtrend2 = zeros(nvars, nobs+1);
    TT2 = zeros(nlag*nvars,1);
    if const>2;
        TT2(1:nvars,:) = AMAT(:,3);
        for i = 2:nobs+1
            HDtrend2_big(:,i) = TT2*((i-1)^2) + Acomp*HDtrend2_big(:,i-1);
            HDtrend2(:,i) = Icomp * HDtrend2_big(:,i);
        end
    end

% Exogenous
    HDexo_big = zeros(nlag*nvars,nobs+1);
    HDexo = zeros(nvars,nobs+1);
    EXO = zeros(nlag*nvars,nvars_ex*(nlag_ex+1));
    if nvars_ex>0;
        VARexo = StructMod.X_EX;
        EXO(1:nvars,:) = AMAT(:,nvars*nlag+const+1:end); % this is c in my notes
        for i = 2:nobs+1
            HDexo_big(:,i) = EXO*VARexo(i-1,:)' + Acomp*HDexo_big(:,i-1);
            HDexo(:,i) = Icomp * HDexo_big(:,i);
        end
    end

% All decompositions must add up to the original data
HDendo = HDinit + HDconst + HDtrend + HDtrend2 + HDexo + sum(HDshock,3);
    
    
    
%% Save and reshape all HDs
%===============================================
HD.shock = zeros(nobs+nlag,nvars,nvars);  % [nobs x shock x var]
    for i=1:nvars
        for j=1:nvars
            HD.shock(:,j,i) = [nan(nlag,1); HDshock(i,2:end,j)'];
        end
    end
HD.init   = [nan(nlag-1,nvars); HDinit(:,1:end)'];    % [nobs x var]
HD.const  = [nan(nlag,nvars);   HDconst(:,2:end)'];   % [nobs x var]
HD.trend  = [nan(nlag,nvars);   HDtrend(:,2:end)'];   % [nobs x var]
HD.trend2 = [nan(nlag,nvars);   HDtrend2(:,2:end)'];  % [nobs x var]
HD.exo    = [nan(nlag,nvars);   HDexo(:,2:end)'];     % [nobs x var]
HD.endo   = [nan(nlag,nvars);   HDendo(:,2:end)'];    % [nobs x var]


%% Plot all HDs
vnames = opt.vnames;

filename = ['./figures/', opt.dataset '_HD_'];
% Initialize HD matrix
[nsteps, nvars, nshocks] = size(HD.shock);


%% Plot
%===============================================
figure('units','normalized','outerposition',[0 0.1 1 0.9])

for ii=1:nvars                
    colormap(parula)
    BarPlot(HD.shock(:,1:nshocks,ii));
    xlim([1 nsteps]);
    title([vnames{ii}], 'FontWeight','bold','FontSize',10); 
    % Save
    FigName = [filename num2str(ii)];
    legend(vnames)
    print('-dpng','-r100',FigName);
    print('-dpdf','-r100',FigName);    
    clf('reset');
end

close all
