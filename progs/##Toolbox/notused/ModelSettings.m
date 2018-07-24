function [ENDO,EXOG,opt] = ModelSettings(opt)
% =======================================================================
% Load data and model settings
% =======================================================================
% LoadDataAndModelSettings(opt)
% -----------------------------------------------------------------------
% INPUT   
%   - opt  : structure with options, see Load Options section below
% -----------------------------------------------------------------------
% OUTPUT
%   - ENDO : endogenous variables, matrix (periods x number of endogenous variables)
%   - EXOG : exogenous variables, matrix (periods x number of exogenous variables)
%   - opt  : structure with options, the following is added to opt:
%            - const: deterministic trend: scalar
%                   - 0 no constant
%                   - 1 constant
%                   - 2 constant and trend
%                   - 3 constant, trend and trend^2
%            - nvars: number of endogenous variables: scalar
%            - nvars_ex: number of exogenous variables: scalar
%            - nlag_ex: number of lags for exogenous variables: scalar
%            - ENDO: matrix of endogenous variables: [periods x nvars]
%            - EXOG: matrix of exogenous variables: [periods x nvars_ex]
%            - dates: cell of strings for dates: [periods x 1]
%            - vnames: cell of strings for names of variables: [1 x nvars]
%            - Restr: structure of restrictions to identify structural shocks.
%                - Rshort:    Short-run restrictions on impact matrix B_0^{-1},
%                             matrix [nvars x nvars]
%                - RlongVAR:  Long-run restrictions in the traditional sense � la Blanchard/Quah,
%                             i.e. on A(1) = inv(eye(nvars)-A1hat-A2hat-...-Aphat)),
%                             matrix [nvars x nvars]
%                - RlongVECM: Long-run restrictions in the VECM model, i.e. on 
%                             Upsylon = (betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)))*B0inv,  
%                             matrix [nvars x nvars]
%                - Please use the following conventions 
%                  - nan: unrestricted
%                  - 0 or any number: zero restriction to specified number
%                  - +1 for positive or -1 for negative sign restriction
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

% Load options
dataset = opt.dataset; % dataset name: string

switch dataset
    %% Canadian
    case 'canadian'
        const = 2;
        [xlsdata, xlstext] = xlsread('datasets/canadian.xlsx','canadian');
        ENDO = xlsdata;
        EXOG=[];
        dates = xlstext(2:end,1);
        dates_freq ='Q';
        vnames = xlstext(1,2:end);
        vnames_ex = [];        
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);        
        nlag_ex = [];        
        Rshort =     [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan 0   nan nan;
                     ];
        RlongVAR  =  [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
        RlongVECM  = [nan 0   0   0;
                      nan nan nan 0;
                      nan nan nan 0;
                      nan nan nan 0;
                     ];
        Rsign  =     [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
    %% Stock and Watson (2001, JEP) data with short-run restrictions
    case 'StockWatson2001'
        const = 2;
        [xlsdata, xlstext] = xlsread('datasets/SW2001_Data.xlsx','Sheet1');
        ENDO = xlsdata;
        EXOG=[];
        dates = xlstext(2:end,1);
        dates_freq ='Q';
        vnames = xlstext(1,2:end);
        vnames_ex = [];
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);
        nlag_ex = [];
        Rshort =     [nan 0   0;
                      nan nan 0;
                      nan nan nan;
                     ];
        RlongVAR  =  [nan nan nan;
                      nan nan nan;
                      nan nan nan;                      
                     ];
        RlongVECM  = [nan nan nan;
                      nan nan nan;
                      nan nan nan;
                     ];
        Rsign  =     [nan nan nan;
                      nan nan nan;
                      nan nan nan;
                     ];
    case 'BlanchardQuah1989_VAR'
        const = 2;
        [xlsdata, xlstext] = xlsread('datasets/BQ1989_Data.xlsx','Sheet1');
        ENDO = xlsdata;
        EXOG=[];
        dates = xlstext(2:end,1);
        dates_freq ='Q';
        vnames = xlstext(1,2:end);
        vnames_ex = [];
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);
        nlag_ex = [];
        Rshort =     [nan nan;
                      nan nan;
                     ];
        RlongVAR  =  [nan 0;
                      nan nan;
                     ];
        RlongVECM  = [nan nan;
                      nan nan;                     
                     ];
        Rsign  =     [nan nan;
                      nan nan;
                     ];
    case 'OilModel1'
        % Quarterly data from FRED
        % GDP: Chain-type Price Index: Index 2009
        % Real Gross Domestic Product: Billions of Chained 2009 Dollars
        % 1973.I-2013.II
        const = 1;
        load datasets/gdpdeflator2.txt; infl=dif(log(gdpdeflator2(56:end,3)))*100;
        load datasets/realgdp2.txt; drgdp=dif(log(realgdp2(56:end,3)))*100;
        % Monthly WTI spot price from Economagic, 1973.1-2007.12
        load datasets/poil.txt; poilm=log(poil(:,3));
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
        dates = dates(1:end-2);
        dates_freq ='Q';
        vnames = {'d(Real Price of Oil)', 'd(GDP Deflator)', 'd(Real GDP)'};
        vnames_ex = [];
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);
        nlag_ex = [];
        Rshort =     [nan 0   0;
                      nan nan 0;
                      nan nan nan;
                     ];
        RlongVAR  =  [nan nan nan;
                      nan nan nan;
                      nan nan nan;
                     ];
        RlongVECM  = [nan nan nan;
                      nan nan nan;
                      nan nan nan;
                     ];
        Rsign  =     [nan nan nan;
                      nan nan nan;
                      nan nan nan;
                     ];
    case 'Keating1992'
        const = 1;
        % Quarterly data from FRED
        % Gross National Product: Chain-type Price Index: Index 2005=100: SA, Q
        % Real Gross Domestic Product: Billions of Chained 2005 Dollars: SAAR, Q
        % Effective federal funds rate, M
        % 1954.IV-2007.IV (exclude unconventional monetary policy after 2007.Q4)
        % Sample reduced to 1959.II-2007.IV
        load datasets/gnpdeflator.txt; infl=dif(log(gnpdeflator(:,3)))*100;
        defl=log(gnpdeflator(:,3));
        load datasets/realgnp.txt; drgdp=dif(log(realgnp(:,3)))*100;
        rgdp=log(realgnp(:,3));
        load datasets/fedfunds.txt;
        irate=[];
        for i=1:3:length(fedfunds(:,3))
          irate=[irate; mean(fedfunds(i:i+2,3))];
        end;    
        % 1959.I-2007.IV reduced to 1959.II-2007.IV
        load datasets/m1.txt;
        m1=log(m1(3:end,3));
        dmoney=[];
        for i=4:3:length(m1)
            dmoney=[dmoney; (m1(i,1)-m1(i-3,1))*100];
        end;
        money=[];
        for i=1:3:length(m1)
            money=[money; m1(i,1)];
        end;
        % Data for Keating (1992) VAR model
        ENDO=[infl(19:end,1) drgdp(19:end,1) irate(19:end,1) dmoney];
        EXOG=[];
                i=1;
        for yyyy = 1959:2007
            for qq = 1:4
                dates{i} = sprintf('%dQ%d',yyyy,qq);
                i=i+1;
            end
        end
        dates = dates(2:end);
        dates_freq ='Q';
        vnames = {'dp','dgnp','i','dm'};
        vnames_ex = [];        
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);        
        nlag_ex = []; 
%         syms b21 b41 b23 b43 b24 b34 b11 b22 b33 b44
%         B0 = [b11 sym(0) sym(0) sym(0);
%              b21  b22 b23    b24;
%              sym(0) sym(0) b33 b34;
%              b41    b41    b43    b44;
%         ];
%         inv(B0) %yields that the first row is restricted to nan 0 0 0
        Rshort =     [nan   0   0   0;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
        RlongVAR  =  [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
        RlongVECM  = [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
        Rsign  =     [nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                      nan nan nan nan;
                     ];
    case 'Gali1999'
        const = 0;
        % Population VAR (1947.I-1998.III) based on Gali (1999), AER, p. 261.
        load datasets/lpmhu.txt; lpmhu=lpmhu(:,2);
        for i=1:length(lpmhu)/3
           hours(i,1)=mean(lpmhu(3*i-2:1:3*i,1));
        end;
        load datasets/gdpq.txt; gdpq=gdpq(:,2);
        ENDO=dif([log(gdpq)-log(hours) log(hours)])*100;
        EXOG=[];
        i=1;
        for yyyy = 1947:1998
            for qq = 1:4
                dates{i} = sprintf('%dQ%d',yyyy,qq);
                i=i+1;
            end
        end
        dates = dates(1:end-2);
        dates_freq ='Q';
        vnames = {'prod','h'};
        vnames_ex = [];
        nvars = size(ENDO,2);       
        nvars_ex = size(EXOG,2);
        nlag_ex = [];
        Rshort =     [nan nan;
                      nan nan;
                     ];
        RlongVAR  =  [nan 0;
                      nan nan;
                     ];
        RlongVECM  = [nan 0;
                      nan nan;                     
                     ];
        Rsign  =     [nan nan;
                      nan nan;
                     ];
end

% Plot dataset
figure;
for i = 1:nvars    
    subplot(nvars,1,i);
    if strcmp(dates_freq,'m')
        plot(datenum(dates,'YYYYmm'),ENDO(:,i));    
        datetick('x','YYYYmm');
    elseif strcmp(dates_freq,'Q')
        plot(datenum(dates,'YYYYQQ'),ENDO(:,i));    
        datetick('x','YYYYQQ');
    end
    title(vnames(i));
end

% Save options
opt.const = const;
opt.dates = dates;
opt.vnames = vnames;
opt.vnames_ex = vnames_ex;
opt.nvars = nvars;
opt.nvars_ex = nvars_ex;
opt.nlag_ex = nlag_ex;
opt.Restr.Rshort = Rshort;
opt.Restr.RlongVAR = RlongVAR;
opt.Restr.RlongVECM = RlongVECM;
opt.Restr.Rsign = Rsign;
