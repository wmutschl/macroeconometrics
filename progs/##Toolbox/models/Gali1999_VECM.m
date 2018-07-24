function [ENDO,EXOG,opt] = Gali1999_VECM(opt)
% =======================================================================
% Replication of Gali (1999) as in Kilian and Lütkepohl (2017, Ch.11)
% The observables are the log of productivity (denoted by prodt) and the 
% log of total employee hours in nonagricultural establishments, 
% averaged to quarterly frequency (denoted by ht). The productivity variable 
% is constructed by subtracting the log of hours from the log of real GDP. 
% Let yt = (prodt, ht) and suppose that both variables are I(1), but
% not cointegrated. 
% Population VAR (1947.I-1998.III) based on Gali (1999), AER, p. 261.
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

load ./datasets/lpmhu.txt;
load ./datasets/gdpq.txt; 
lpmhu=lpmhu(:,2);
for i=1:length(lpmhu)/3
   hours(i,1)=mean(lpmhu(3*i-2:1:3*i,1));
end;
gdpq=gdpq(:,2);

ENDO=([log(gdpq)-log(hours) log(hours)])*100;
EXOG=[];
i=1;
for yyyy = 1947:1998
    for qq = 1:4
        dates{i} = sprintf('%dQ%d',yyyy,qq);
        i=i+1;
    end
end
opt.dates = dates(1:end-1);
opt.dates_freq ='Q';
opt.vnames = {'prod','h'};
opt.vnames_ex = [];


opt.Rshort = [nan nan;
              nan nan;
             ];
opt.Rshorth = {};
% Only technology shock has a long-run effect on the level of productivity
opt.Rlong  = [nan 0;
              nan nan;
             ];
opt.Rsign  = [nan nan;
              nan nan;
             ];
opt.Rsignh = {};

%% Model settings
opt.const = 1;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'ML';   % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H1';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 0;         % Scalar, values are: 1 plot cumulated IRFs, 0 else
