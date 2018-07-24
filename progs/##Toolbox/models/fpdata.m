function [ENDO,EXOG,opt] = fpdata(opt)
% =======================================================================
% Replication of: 

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
%           - Rsignh:  Cell of Matrix of sign restrictions, [nvars x nvars]

% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu
x = csvread('./datasets/fpdata5.csv', 1, 1);


%--------------------------------------------------------data construction 

%construct nonfarm business sector inflation rate, losing one observation
p=x(:,8);
infl=4*(log(p(2:end))-log(p(1:end-1)));
x=x(2:end,:); % now everything is 1948.2 - 2013.4

ramey=x(:,2);
tbr=x(:,3);
barrotax=x(:,4);
hours=x(:,5);
govsaving=x(:,6);
y=x(:,7);
p=x(:,8);
dtfp=x(:,9);
g=x(:,10);
gdpr=x(:,11);
gdpp=x(:,12);
invrtotal=x(:,13); %total private investment, includes housing
dutil=x(:,14);
govcons=x(:,15);
govint=x(:,16);
%
tax=(x(:,17)+x(:,19)-x(:,18))./(gdpp/100); %real net taxes, as described in text
%
%nonresidential investment:
invr=x(:,21); %real private nonresidential investment, quantity index
invn=x(:,22); %nominal private nonresidential investment
    %observation 245 is 2009q1, convert quantitiy index to real value in
    %2009 prices
invr=(invr./100).*invn(245);

%interest rate:
tbr=tbr/100;

%the government deficit, fraction of nominal gdp:
gdpn=gdpr.*gdpp/100;
def=-govsaving./gdpn;

    % dtfp and dutil are given by Fernald as annualized growth rates, here I construct the levels 
    % of TFP and UTIL from this information to be consistent with the
    % levels formulation of the rest of the variables 
    % I iterate on dtfp_t=400*(log tfp_t - log tpf_t-1) using 1 as starting
    % value
    tfp=zeros(length(dtfp),1);
    util=zeros(length(dutil),1);
        tfp(1)=exp(dtfp(1)/400);
        util(1)=exp(dutil(1)/400);
    for t=2:length(tfp)
        tfp(t)=exp(dtfp(t)/400+log(tfp(t-1)));
        util(t)=exp(dutil(t)/400+log(util(t-1)));
    end
    
spec = 1;
% Spec 1-3:   Recursive Perotti identification scheme, direct labprod (log(y/h)).
% Spec 4-6:   Recursive Ramey identification scheme, direct labprod (log(y/h)).
% Spec 7-9:   Recursive Perotti identification scheme, indirect labprod (hours).
% Spec 10-12: Recursive Ramey identification scheme, indirect lkabprod (hours).
% Spec 13-15: Sign Restrictions Specification

if spec == 1
    ENDO = [ log(g), log(y), log(y./hours), tbr, log(tax), infl, def ];
    opt.vnames = {'Govt spending', 'GDP', 'Labprod', 'T-bill rate', 'Barro tax', 'Inflation', 'Deficit'};
elseif spec == 2
    ENDO =  [log(g), log(y), log(y./hours), tbr, log(tax), infl, def, log(util), log(invr)];
    opt.vnames =  {'Govt spending', 'GDP', 'Labprod', 'T-bill rate', 'Barro tax', 'Inflation', 'Deficit', 'Util rate', 'Investment'};
elseif spec == 3
    ENDO = [log(g), log(y), log(y./hours), tbr, log(tax), infl, def, log(invr)];
    opt.vnames =  {'Govt spending', 'GDP', 'Labprod', 'T-bill rate', 'Barro tax', 'Inflation', 'Deficit', 'Investment'};

elseif spec == 4
    data = [ramey, log(g), log(y), log(y./hours), tbr, log(tax), infl, def ];
elseif spec == 5
    data =  [ramey, log(g), log(y), log(y./hours), tbr, log(tax), infl, def, log(util), log(invr)];
elseif spec == 6
    data =  [ramey, log(g), log(y), log(y./hours), tbr, log(tax), infl, def, log(util), log(invr)];

elseif spec == 7
    data = [log(g), log(y), tbr, log(tax), log(hours), infl, def ];
elseif spec == 8
    data =  [log(g), log(y), tbr, log(tax), log(hours), infl, def, log(util), log(invr)];
elseif spec == 9
    data =  [log(g), log(y), tbr, log(tax), log(hours), infl, def, log(invr)];

elseif spec == 10
    data = [ramey, log(g), log(y), tbr, log(tax), log(hours), infl, def ];
elseif spec == 11
    data =  [ramey, log(g), log(y), tbr, log(tax), log(hours), infl, def, log(util), log(invr)];
elseif spec == 12
    data =  [ramey, log(g), log(y), tbr, log(tax), log(hours), infl, def, log(invr)];

elseif spec == 13
    %- Without util:
    ENDO=[log(g) log(y) tbr log(tax) log(hours) infl def];
    opt.vnames = {'g' 'y' 'tbr' 'tax' 'hours' 'infl' 'def'};
elseif spec == 14
    %- With util:
    ENDO=[log(g) log(y) tbr log(tax) log(hours) infl def log(util) log(invr)];
    opt.vnames = {'g' 'y' 'tbr' 'tax' 'hours' 'infl' 'def' 'util' 'invr'};
elseif spec == 15
    %- With tfp:
    ENDO=[log(g) log(y) tbr barrotax log(hours) infl def log(util) log(invr) log(tfp)];
    opt.vnames = {'g' 'y' 'tbr' 'tax' 'hours' 'infl' 'def' 'util' 'invr' 'tfp'};
end

EXOG=[];
opt.dates = [];
opt.dates_freq = 'Q';
opt.vnames_ex = [];

opt.Rshort    = [nan nan nan nan nan;
                 nan nan nan nan nan;
                 nan nan nan nan nan;
                 nan nan nan nan nan;
                 nan nan nan nan nan;
                 ];
opt.Rshorth = {}; % if there are no short-run restrictions at horizon h

opt.Rlong  = [nan nan nan nan nan;
              0   0   nan 0   0;
              nan nan nan nan nan;
              nan nan nan nan nan;
              nan nan nan nan nan;
             ];

% Sign restrictions on impact
opt.Rsign  = [+1  +1  -1  -1 nan;
              -1  +1  +1  +1 nan;
              nan nan nan +1 nan;
              -1  +1  -1  -1 nan;
              nan nan nan -1 nan;
             ];
opt.Rsignh = {};

%% Model settings
opt.const = 3;             % Set deterministic trend: 0 no constant; 1 constant; 2 constant and trend; 3 constant, trend and trend^2
opt.EstMethodVAR = 'OLS';  % Estimation method for VAR models: 'OLS': ordinary least squares, 'ML': maximum likelihood
opt.JOHmodel = 'H1';       % String specifying the form of the deterministic components of the VEC(nlag-1) model, see jcitest.m. Values are
                           % -'H2':    A*B'*y(t-1). There are no intercepts or trends in the cointegrating relations and there are no trends in the data.
                           % -'H1*':   A*(B'*y(t-1) + c0). There are intercepts in the cointegrating relations and there are no trends in the data.
                           % -'H1':    A*(B'*y(t-1) + c0) + c1. There are intercepts in the cointegrating relations and there are linear trends in the data.
                           % -'H*':    A*(B'*y(t-1) + c0 + d0*t) + c1. There are intercepts and linear trends in the cointegrating relations and there are linear trends in the data.
                           % -'H':     A*(B'*y(t-1) + c0 + d0*t) + c1 + d1*t. There are intercepts and linear trends in the cointegrating relations and there are quadratic trends in the data.
opt.IRFcumsum = 0;         % Scalar, values are: 1 plot cumulated IRFs, 0 else
