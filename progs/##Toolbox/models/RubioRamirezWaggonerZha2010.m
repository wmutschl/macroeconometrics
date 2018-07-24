function [ENDO,EXOG,opt] = RubioRamirezWaggonerZha2010(opt)
% =======================================================================
% Replication of Rubio-Ramirez, Waggoner and Zha (2010) as in Kilian and 
% Lütkepohl (2017, Ch.11).
% Let zt = (d(gnpt), it,d(pt)) ? I(0), where gnpt denotes the log of
% U.S. real GNP, pt the corresponding GNP deflator in logs, and it the federal
% funds rate, averaged by quarter. The estimation period is restricted to 1954q4-
% 2007q4 in order to exclude the recent period of unconventional monetary policy measures.
% Quarterly data from FRED
% Gross National Product: Chain-type Price Index: Index 2005=100: SA, Q
% Real Gross Domestic Product: Billions of Chained 2005 Dollars: SAAR, Q
% Effective federal funds rate, M
% 1954.IV-2007.IV (exclude unconventional monetary policy after 2007.Q4)
% =======================================================================
% Please use the following conventions for the restrictions
% - nan: unrestricted
% - 0 or any number: zero restriction to specified number (any number does not work for RWZ)
% - +1 for positive or -1 for negative sign restriction
% =======================================================================
% OUTPUT
%   - ENDO : endogenous variables, matrix (periods x number of endogenous variables)
%   - EXOG : exogenous variables, matrix (periods x number of exogenous variables)
%   - opt  : structure with options:
%            - dates: cell of strings for dates: [periods x 1]
%            - dates_freq: 'Q' for quarterly data, 'm' for monthly data
%            - vnames: cell of strings for names of endogenous variables: [1 x nvars]
%            - vnames_ex: cell of strings for names of exogenous variables: [1 x nvars_ex]
%            - Rshort: Short-run restrictions on impact matrix B_0^{-1}, matrix [nvars x nvars]
%            - Rlong:  Long-run restrictions, matrix [nvars x nvars]
%                      - for VAR on A(1)*B0inv = inv(eye(nvars)-A1hat-A2hat-...-Aphat))B0inv
%                      - for VECM on Upsylon*B0inv = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))*B0inv,  
%           - Rsign:  Sign restrictions, matrix [nvars x nvars]
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

load ./datasets/gnpdeflator.txt; 
load ./datasets/realgnp.txt; 
load ./datasets/fedfunds.txt;
infl=dif(log(gnpdeflator(:,3)))*100;
drgdp=dif(log(realgnp(:,3)))*100;
irate=[];
for i=1:3:length(fedfunds(:,3))
  irate=[irate; mean(fedfunds(i:i+2,3))];
end;    
ENDO=[drgdp irate infl];
EXOG=[];
i=1;
for yyyy = 1954:2007
    for qq = 1:4
        dates{i} = sprintf('%dQ%d',yyyy,qq);
        i=i+1;
    end
end
opt.dates = dates(4:end);
opt.dates_freq ='Q';
opt.vnames = {'d(gnpt)', 'it', 'd(pt)'};
opt.vnames_ex = [];

% Monetary policy shock does not affect real GNP in the current quarter
opt.Rshort = [0   nan nan;
              nan nan nan;
              nan nan nan;
             ];
opt.Rshorth = {};
% Only the aggregate supply shock affects real GNP in the long-run
opt.Rlong  = [0   0   nan;
              nan nan nan;
              nan nan nan;
             ];
opt.Rsign  = [nan nan nan;
              nan nan nan;
              nan nan nan;
             ];
opt.Rsignh = {};

%% Model settings
opt.const = 1;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'OLS';  % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H*';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 0;         % Scalar, values are: 1 plot cumulated IRFs, 0 else
