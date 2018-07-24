function result=chapter4_unitRootTesting()
%% PURPOSE: Test data for unit root
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The code is not written for computational efficiency or
% elegance. The book: 
%
% "Applied Time Series for Macroeconomics"
% Gyldendal Akademisk 2014
% by Hilde C. Bjørnland and Leif A. Thorsrud 
%
% provides details. Please refer to the book if the code(s) are used for 
% research of commercial purposes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Settings

% Without trend and constant p=-1, with constant p=0, with constant and trend p=1
p=0:1; 
% Number of lags to use in ADF test
nlag=[2 4 6];

%% Load data
%names={'Oslo stock exchange','3-month interest rate','Trade weighted exchange rate','CPI','GDP','Unemployment'};
sname={'osebx','inter3m','i44','cpi','gdp','unemp'};

datam=xlsread('./timeSeriesData.xlsx','m')';
dataq=xlsread('./timeSeriesData.xlsx','q')';
data=[mat2cell(datam(:,2:end),size(datam,1),ones(1,3)) mat2cell(dataq(:,2:end),size(dataq,1),ones(1,3))];

%% Do unit root testing employing ADF
nump=numel(p);
numnlag=numel(nlag);

% result will be a structure containing in ith element test results for
% different trends constants etc., in the jth element test results for
% different lag lenghts and finally in the h element the data. 
for h=1:numel(data)
    for i=1:nump
        for j=1:numnlag
            d=data{1,h};
            d=d(~isnan(d));            
            result(i,j,h)=unitRootTest(d,p(i),nlag(j),'vname',sname{h});
        end;
    end;
end;






