function [ENDO,EXOG,opt] = KingRebeloPlosserWatson1991_VAR(opt)
% =======================================================================
% Replication of: King, Rebelo, Plosser and Watson (1991,AER).
% Baseline Three-Variable Model: y = [d(gnp) (c-gnp) (inv-gnp)]'
% gnp: US real output
% c  : US real consumption
% inv: US real investment
% Productivity shock with stochastic trend affects all three variables
% Balanced growth implies that (c-gnp)~I(0) and (inv-gnp)~I(0)
% Estimation of Vector Autoregressive Model
% Purpose is to study the effect of the productivity shock, other two shocks are identified arbitrarily
% Important: Need to cumulate IRFs!
% =======================================================================
% Please use the following conventions for the restrictions
% - nan: unrestricted
% - 0 or any number: zero restriction to specified number (any number does not work for RWZ algorithm)
% - +1 for positive or -1 for negative sign restriction
% =======================================================================
% INPUTS
%   - opt : structure with options to  be updated
% =======================================================================
% OUTPUTS
%   - ENDO : endogenous variables, matrix (periods x number of endogenous variables)
%   - EXOG : exogenous variables, matrix (periods x number of exogenous variables)
%   - opt  : structure with options:
%            - dates: cell of strings for dates, [periods x 1]
%            - dates_freq: 'Q' for quarterly data, 'm' for monthly data
%            - vnames: cell of strings for names of endogenous variables, [1 x nvars]
%            - vnames_ex: cell of strings for names of exogenous variables, [1 x nvars_ex]
%            - Rshort: Matrix of short-run restrictions on impact matrix B_0^{-1}, [nvars x nvars]
%            - Rlong:  Matrix of long-run restrictions, [nvars x nvars]
%                      - for VAR on A(1)*B0inv = inv(eye(nvars)-A1hat-A2hat-...-Aphat))B0inv
%                      - for VECM on Upsylon*B0inv = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))*B0inv,  
%           - Rsign:   Matrix of sign restrictions, [nvars x nvars]
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu%/* ***** Dates for Data set ****** */
y=load('./datasets/KPSW_ciy.txt');  % load data from file
gdp = y(:,3);
c = y(:,1);
inv = y(:,2);
dgdp = gdp(2:end,:) - gdp(1:end-1,:);
c_gnp = c(2:end,:) - gdp(2:end,:);
inv_gnp = inv(2:end,:) - gdp(2:end,:);

ENDO=[dgdp c_gnp inv_gnp];
EXOG=[];
i=1;
for yyyy = 1947:1988
    for qq = 1:4
        dates{i} = sprintf('%dQ%d',yyyy,qq);
        i=i+1;
    end
end
opt.dates = dates(2:end);
opt.dates_freq ='Q';
opt.vnames = {'dgnp', 'c-gdp', 'inv-gdp'};
opt.vnames_ex = [];

opt.Rshort  = [nan nan nan;
               nan nan nan;
               nan nan nan;
              ];
opt.Rshorth = {}; % if there are no short-run restrictions at horizon h

% only productivity shocks have long-run effect on output, non-productivity shocks
% are arbitrarily identified, e.g. lower cholesky
opt.Rlong  = [nan 0   0;
              nan nan 0;
              nan nan nan;
             ];

% Sign restrictions on impact
opt.Rsign  = [nan nan nan;
              nan nan nan;
              nan nan nan;
             ];
% Sign restrictions at horizon h
opt.Rsignh  = {};

%% Model settings
opt.const = 1;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'OLS';  % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H1';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 1;         % Scalar, values are: 1 plot cumulated IRFs, 0 else
