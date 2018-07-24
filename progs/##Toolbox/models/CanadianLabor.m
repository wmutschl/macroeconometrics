function [ENDO,EXOG,opt] = CanadianLabor(opt)
% =======================================================================
% Replication of: Breitung, Brüggemann, Lütkepohl (2004) - Structural Vector 
% Autoregressive Modeling and Impulse Responses.
% Structural VECM for Canadian labor market data specified for labor
% productivity, employment, unemployment rate, and real wages. 
% Data is quarterly, seasonally adjusted constructed from data
% obtained from the OECD database for the period from 1980Q1 to 2000Q4
% Productivity (gdp - e) is constructed by subtracting the log of employment (e)
% from the log real GDP (gdp), U denotes the unemployment rate, and w - p is
% the log of a real wage index. 
% The time series vector is y t = [(gdp_t - e_t), e_t , U_t , (w_t - p_t)].
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
%            - Rshorth: Cell of matrix of short-run restrictions on impact matrix at horizon h, i.e. on (J*Acomp^h*J')*B_0^{-1}, cell of [nvars x nvars] matrix
%                       For horizon h specify Rshorth{h} = [nvars x nvars] matrix; if there are no short-run restrictions at horizon h, please set Rshorth = {}, 
%            - Rlong:  Matrix of long-run restrictions, [nvars x nvars]
%                      - for VAR on A(1)*B0inv = inv(eye(nvars)-A1hat-A2hat-...-Aphat))B0inv
%                      - for VECM on Upsylon*B0inv = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))*B0inv,  
%           - Rsign:   Matrix of sign restrictions, [nvars x nvars]
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

[xlsdata, xlstext] = xlsread('./datasets/canadian.xlsx','canadian');
ENDO = xlsdata;
EXOG=[];
opt.dates = xlstext(2:end,1);
opt.dates_freq ='Q';
opt.vnames = xlstext(1,2:end);
opt.vnames_ex = [];
% Labor demand shocks do not affect real wages on impact, zero in (4,2)
opt.Rshort    = [nan nan nan nan;
                 nan nan nan nan;
                 nan nan nan nan;
                 nan 0   nan nan;
                 ];

% Short-run restrictions at horizon h
opt.Rshorth = {}; % if there are no short-run restrictions at horizon h
%h = 2;
%opt.Rshorth{h} = [nan nan nan nan;
%                 nan nan nan nan;
%                 nan nan nan nan;
%                 nan 0   nan nan;
%                 ];
% Wage shocks are transitory, i.e. have no long-run impact (last column of zeros)
% Due to constant returns to scale, productivity is only driven by technology shocks, zeros in (1,2), (1,3) and (1,4)
opt.Rlong  = [nan 0   0   0;
              nan nan nan 0;
              nan nan nan 0;
              nan nan nan 0;
             ];
% Sign restrictions on impact
opt.Rsign  = [nan nan nan nan;
              nan nan nan nan;
              nan nan nan nan;
              nan nan nan nan;
             ];
% Sign restrictions at horizon h
opt.Rsignh  = {};

%% Model settings
opt.const = 2;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'OLS';  % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H*';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 0;         % Scalar, values are: 1 plot cumulated IRFs, 0 else
