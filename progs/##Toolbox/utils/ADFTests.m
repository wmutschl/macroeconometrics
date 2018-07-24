function ADFTests(ENDO,opt)
% =======================================================================
% Perform and display summary of Augmented Dickey Fuller tests for variables
% in levels as well as in first differences
% =======================================================================
% ADFTests(ENDO,opt)
% -----------------------------------------------------------------------
% INPUTS
%   - ENDO : matrix [periods x number of variables]
%   - opt  : structure with options, see load options section below
% -----------------------------------------------------------------------
% CALLS
%   - adftest.m: Augmented Dickey-Fuller test for a unit root
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu

% Load Options
nvars = size(ENDO,2);           % number of variables: scalar
vnames = opt.vnames;            % names of variables: (nvars x 1) cell of strings
maxlags = opt.ADFmaxlags;       % maximum number of lags to consider: scalar
optlagcrit = opt.ADFoptlagcrit; % criterium to choose optimal lag: string, values: 'AIC','BIC','HQC'

for i = 1:nvars
    for j=1:2 
        if j == 1    % Perform tests for variables in levels
            Y = ENDO(:,i);
            nameY = vnames{i};
        elseif j== 2 % Perform tests for variables in first differences
            Y = diff(ENDO(:,i));
            nameY = strcat('d(',vnames{i},')');
        end
        % Basic autoregressive model, see adftest.m
        [hAR,pValueAR,statAR,cValueAR,regAR] = adftest(Y,'model','AR','lags',0:maxlags);
        [~,iAR]=min([regAR.(optlagcrit)]); % get optimal lag according to criterium
        % Autoregressive model with drift, see adftest.m
        [hARD,pValueARD,statARD,cValueARD,regARD] = adftest(Y,'model','ARD','lags',0:maxlags);
        [~,iARD]=min([regARD.(optlagcrit)]); % get optimal lag according to criterium
        % Trend stationary autoregressive model, see adftest.m
        [hTS,pValueTS,statTS,cValueTS,regTS] = adftest(Y,'model','TS','lags',0:maxlags); 
        [~,iTS]=min([regTS.(optlagcrit)]); % % get optimal lag according to criterium    
        % Store results, i.e. optimal lages, test statistics and p values
        results=[iAR-1 iARD-1 iTS-1;
                 statAR(iAR), statARD(iARD), statTS(iTS);
                 pValueAR(iAR), pValueARD(iARD), pValueTS(iTS);
                ];
        % Display summary of test results
        fprintf('*********************************\n');
        fprintf('*** ADF TEST RESULTS FOR %s:\n',nameY);
        fprintf('*********************************\n');
        disp(array2table(results,'rowNames',{sprintf('Optimal Lags (%s)',optlagcrit),'Test Value','p Value'},'VariableNames',{'AR','ARD','TS'}))
        if isempty(find(hAR==1)) == 0;  fprintf('   Reject Null for %s in AR model for lags: %s\n',nameY,mat2str(find(hAR==1)-1));   end
        if isempty(find(hARD==1)) == 0; fprintf('   Reject Null for %s in ARD model for lags: %s\n',nameY,mat2str(find(hARD==1)-1)); end
        if isempty(find(hTS==1)) == 0;  fprintf('   Reject Null for %s in TS model for lags: %s\n',nameY,mat2str(find(hTS==1)-1));   end
        fprintf('\n\n');
    end
end