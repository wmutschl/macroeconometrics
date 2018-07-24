function [ENDO,EXOG,opt] = KilianLuetkepohl2017USOil(opt)
% =======================================================================
% Replication of the structural VAR for the US oil in Kilian and Luetkepohl (2017, Ch. 9)
% Quarterly data from FRED
% GDP: Chain-type Price Index: Index 2009, 1973.I-2013.II
% Real Gross Domestic Product: Billions of Chained 2009 Dollars, 1973.I-2013.II
% Monthly WTI spot price from Economagic, 1973.1-2007.12
% =======================================================================
% Please use the following conventions for the restrictions
% - nan: unrestricted
% - 0 or any number: zero restriction to specified number (any number does not work for RWZ)
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
%           - Rsignh:  Cell of Matrix of sign restrictions, [nvars x nvars]
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

load ./datasets/gdpdeflator2.txt; 
load ./datasets/realgdp2.txt; 
load ./datasets/poil.txt; 
poilm=log(poil(:,3));
drgdp=dif(log(realgdp2(56:end,3)))*100;
infl=dif(log(gdpdeflator2(56:end,3)))*100;
drpoil=(poilm(4:3:end)-poilm(1:3:end-3))*100-infl;

ENDO=[drpoil infl drgdp];        
EXOG=[];
i=1;
for yyyy = 1973:2013
    for qq = 1:4            
        dates{i} = sprintf('%dQ%d',yyyy,qq);
        i=i+1;
    end
end
opt.dates = dates(1:end-2);
opt.dates_freq ='Q';
opt.vnames = {'d(Real Price of Oil)', 'd(GDP Deflator)', 'd(Real GDP)'};
opt.vnames_ex = [];

% Recursive structure: predetermined real oil prices
opt.Rshort =  [nan 0   0;
               nan nan 0;
               nan nan nan;
              ];
opt.Rshorth = {}; % if there are no short-run restrictions at horizon h

opt.Rlong  =  [nan nan nan;
               nan nan nan;
               nan nan nan;
              ];

opt.Rsign  =  [nan nan nan;
               nan nan nan;
               nan nan nan;
              ];
opt.Rsignh  = {};

%% Model settings
opt.const = 1;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'OLS';  % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H*';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 1;         % Scalar, values are: 1 plot cumulated IRFs, 0 else

