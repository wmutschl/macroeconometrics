%% Toolbox for Structural Vector Autogressive (SVAR) and Structural Vector Error Correction (SVECM) Models 
% =======================================================================
% Description:
% This suite of Matlab functions provides a comprehensive toolbox to 
% estimate structural vector autoregressive models (SVAR) as well as 
% structural vector error correction models (SVECM).
% Three identification schemes may be specified (separately or combined): 
% short-run restrictions, long-run restrictions, long-run restrictions due 
% to cointegration and sign restrictions.
% Local and global identification criteria are also checked.
%
% To this end, the toolbox does the following:
% 1) Summary of ADF Tests for variables in levels and first differences
% 2) Cointegration analysis using Johansen's maxeig or trace test
% 3) Estimation of reduced-form model: VAR is estimated via OLS, VECM via 
%    Johansen's reduced-rank Maximum Likelihood
% 4) In case of just-identified or over-identified models:
%      Estimation of the structural model using one of four algorithms:
%      - Cholesky decomposition for recursive models
%      - Numerical optimization using Matlab's fsolve
%      - QR algorithm of Rubio-Ramirez,Waggoner, Zha (2010) or Binning (2013)
%      - Scoring algorithm of Amisano and Giannini (1997)
%      The normalization rule imposes that the sign of the diagonal elements
%      of B0inv are positive
% 5) If the model is underidentified using only exclusion restrictions, 
%    the toolbox uses sign restrictions to narrow down the set of valid 
%    models using the QR algorithm of Binning (2013)
% 6) Impulse Responses and bootstrapped bands are computed and saved to 
%    files
% 7) Forecast Error Variance Decompositions and bootstrapped bands are 
%    computed and saved to files
% 8) Historical Decompositions are computed and saved to files
% -----------------------------------------------------------------------
% USAGE:
% Please set all all model settings and restrictions in separate m files, 
% see the examples in the "model" folder. All other options are set at the 
% beginning of this script file or can be chosen interactively
% -----------------------------------------------------------------------
% NOTES:
% This toolbox encompasses modified versions of codes provided by Andrew Binning, 
% Ambrogio Cesa Bianchi, Lutz Kilian, James P.LeSage, and JMulTi 
% contriubutors (Joerg Breitung, Ralf Brueggemann, Helmut Herwartz, Markus 
% Kraetzig, Helmut Luetkepohl, Timo Teraesvirta, and Rolf Tschernig)
% -----------------------------------------------------------------------
% AUTHOR:
% Willi Mutschler, May 2017
% willi@mutschler.eu
% -----------------------------------------------------------------------
% Things in process or not working yet:
% - Solve Error in PlotData.m
% - make option to change normalization rule
% - Works only with B-model, 
% - Implement function that checks order and rank condition of identification for VAR & VECM
% - Implement Sign Restrictions for exactly identified models
% - Implement Sign Restrictions for underidentified models, (still bootstrap in IRF functions or redundant???)
% =======================================================================

%% House keeping
clearvars; clearvars global; close all; clc;
addpath('utils/','models/','-begin');

%% Set Options
% =======================================================================
% Set datasets: Please see corresponding m files in models folder for details
%   - 'BinningExample': VAR(2)
%   - 'BlanchardQuah1989': VAR(4)
%   - 'CanadianLabor': VECM(2) with cointegration rank r=1
%   - 'Gali1999_VAR': VAR(4)
%   - 'Gali1999_VECM': VECM(4) with cointegratrion rank r=0
%   - 'KilianLuetkepohl2017USOil': VAR(4)
%   - 'KingRebeloPlosserWatson1991_VAR': VAR(2)
%   - 'KingRebeloPlosserWatson1991_VECM': VECM(1) with cointegration rank r = 2
%   - 'RubioRamirezWaggonerZha2010': VAR(4)
%   - 'StockWatson2001': VAR(4)
%   - 'fpdata': VAR(4)

opt.dataset = 'BlanchardQuah1989';

opt.prelimtesting = 1;     % 1 to do preliminary testing, i.e. ADF, lag selection and cointegration tests
opt.algorithm = 'Scoring';     % String for algorithm. Values are:
                           % - 'Cholesky' to use Cholesky decomposition for recursive models
                           % - 'Optimization' to use numerical optimization routine fsolve
                           % - 'RWZ' in the fashion of Binning (2013), Arias, Rubio-Ramirez and Waggoner (2016) or Rubio-Ramirez, Waggoner and Zha (2010)
                           % - 'Scoring' in the fashion of Amisano and Giannini (1997), Breitung, Br�ggemann und L�tkepohl (2004) or Gauss code of JMulTi
opt.StartValueMethod = 1;  % Scalar for starting values of algorithms
                           % - RWZ & Optimization: 1 Cholesky decomposition 2 square root of covariance matrix of reduced-form residuals
                           % - Scoring: 1 draw randomly, 2 fixed values (see StructuralForm.m)
opt.ADFmaxlags = 5;        % Scalar for maximum number of lags considered in ADF Tests
opt.ADFoptlagcrit = 'BIC'; % String for criteria used to select lags in ADF Tests, values are: 'AIC','BIC','HQ'
opt.maxlags = 5;           % Scalar for maximum number of lags to consider for lag selection test
opt.JOHtest = 'trace';     % String indicating the type of cointegration test to be performed. Values are 'trace' or 'maxeig', see jcitest.m for further details
opt.nsteps = 20;           % Scalar for number of steps for IRFs and FEVDs
opt.impact = 0;            % Scalar for size of the shock for IRFs (0 for one standard deviation, 1 for 1)
opt.ndraws = 100;          % Scalar for draws for bootstrap and sign restrictions
opt.pctg   = 95;           % Sclar for confidence bands for bootstrap
opt.method = 'bs';         % String for methodology for error bands, values are 'bs' for sampling with replacement of residuals bootstrap, 'wild' for wild sampling of residuals bootstrap

%% Specification and Reduced-Form estimation
% =======================================================================
[ENDO,EXOG,opt] = feval(opt.dataset,opt); %load data and model settings
if isempty(EXOG) == 0
    opt.nlag_ex = input('Please set the lag order for the exogenous variables: ');
else
    opt.nlag_ex = [];
end
%PlotData(ENDO,EXOG,opt)
if opt.prelimtesting
    ADFTests(ENDO,opt);                                % Perform ADF Unit Root tests for each ENDO variable
    [opt.nlag,VAR] = LagOrderSelection(ENDO,EXOG,opt); % Lag order selection and OLS estimation of VAR
    [opt.coint_r,VECM] = CointTest(ENDO,EXOG,opt);     % Johansen cointegration test and ML estimation of VECM    
    opt.model = input(sprintf('Input 1 for VAR(%d), or 2 for VECM(%d) with cointegration rank r=%d: ',opt.nlag,opt.nlag-1,opt.coint_r));
    if opt.model == 1 
        ReducedForm = VAR;
    elseif opt.model == 2
        ReducedForm = VECM;
    end
else
    opt.nlag = input('Please set the lag order for VAR(nlags) or VECM(nlags-1): ');
    opt.model = input(sprintf('Input 1 for VAR(%d), or 2 for VECM(%d) : ',opt.nlag,opt.nlag-1));
    if opt.model == 2
        opt.coint_r = input(sprintf('Please set the cointegration rank r equal to: '));
    end
    ReducedForm = EstimReducedForm(ENDO,EXOG,opt); % Estimate reduced-form VAR or VECM
end
fprintf('Estimate of covariance-matrix of reduced-form shocks:\n'); disp(ReducedForm.SIGMAUHAT);

%% Structural Estimation
% =======================================================================
[StructMod,opt] = StructuralForm(ReducedForm,opt); % Estimate structural model
%% Impulse Responses
% =======================================================================
if StructMod.SignRestr == 0
    fprintf('Estimate of short-run impact matrix:\n'); disp(StructMod.B0inv_point);
    fprintf('Estimate of long-run impact matrix:\n');  disp(StructMod.LRMat_point);

    % Compute IRF
    IRF.point = IRFs(StructMod,opt);
    IRFplot(opt,IRF.point,[],[]);
    % Compute error bands
    [IRF.INF,IRF.SUP,IRF.MED,IRF.BAR] = IRFband(StructMod,opt);
    % Plot and save to files into folder /figures/
    IRFplot(opt,IRF.MED,IRF.INF,IRF.SUP);


    %% Forecast Error Variance Decomposition
    % =======================================================================
    % Compute FEVD
    FEVD.point = FEVDs(StructMod,opt);
    FEVDplot(opt,FEVD.point,[],[]);
    % Compute error bands
    [FEVD.INF,FEVD.SUP,FEVD.MED] = FEVDband(StructMod,opt);
    % Plot and save to files into folder /figures/
    FEVDplot(opt,FEVD.MED,FEVD.INF,FEVD.SUP);

    %% Historical Decomposition
    % =======================================================================
    % Compute, plot and save figures to files into folder /figures/
    HD = HistDecomp(StructMod,opt);
elseif StructMod.SignRestr == 1
    % 1st dimension = variable
    % 2nd dimension = time
    % 3rd dimension = draw
    % 4th dimension = shock
    plot(squeeze(StructMod.IRFs(1,1:20,:,1)));
end