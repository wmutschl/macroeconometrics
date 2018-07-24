function [coint_r,VECM] = CointTest(ENDO,EXOG,opt)
% =======================================================================
% Perform and display summary of cointegration tests
% =======================================================================
% [coint_r,VECM] = CointTest(ENDO,EXOG,opt)
% -----------------------------------------------------------------------
% INPUTS
%   - ENDO : matrix of endogenous variables, [periods x number of endogenous variables)]
%   - EXOG : matrix of exogenous variables, [periods x number of exogenous variables]. Please note: jcitest does not support exogenous variables
%   - opt  : structure with options to hand over to EstimReducedForm.m
% -----------------------------------------------------------------------
% OUTPUT
%   - coint_r: cointegration rank choosen by the user, scalar
%   - VECM:  structure including VECM(nlag-1) estimation results for cointegration rank coint_r
% -----------------------------------------------------------------------
% CALL
%   - EstimReducedForm.m to estimate reduced-form VECM model
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

if isempty(EXOG) == 0
    error('jcitest does not support exogenous variables, discarding them.')
end
% Call jcitest to perform and display the Johansen cointegration test
opt.model = 2; % Estimate VECM not VAR model
opt.JOHCointTest = 1; % Display results of tests
VECM = EstimReducedForm(ENDO,EXOG,opt);
coint_r = VECM.coint_r; % cointegration rank
% DiagnosticTests(ENDO,EXOG,opt); % Diagnostic Tests against autocorrelation, nonnormality and ARCH effects in VAR residuals not yet