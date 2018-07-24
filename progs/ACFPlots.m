function [r_k] = ACFPlots(y,pmax,alph)
% =======================================================================
% Computes and plots the empirical autocorrelation function
% c_k = 1/T*\sum_{t=k+1}^T (y_t-\bar{y})(y_{t-k}-\bar{y})
% r_k = c_k/c0
% =======================================================================
% ACFPlots(y,pmax,alph)
% -----------------------------------------------------------------------
% INPUTS
%   - y    : Vector of data. [periods x 1]
%   - pmax : Maximum number of lags to plot. [scalar]
%   - alph : Significance level for asymptotic bands. [scalar]
% -----------------------------------------------------------------------
% OUTPUTS
%   - rk   : Sample autocorrelation coefficient. [1 x maxlags]
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================

T=size(y,1);                            % get number of periodes
y_demeaned = bsxfun(@minus,y,mean(y));  % put y in deviations from mean
r_k = nan(1,pmax);                      % initialize output vector

% Compute variance
c0 = 1/T*(y_demeaned' * y_demeaned);
% Compute autocorrelations
for k=1:pmax
    r_k(1,k) = 1/(T*c0) * (y_demeaned(1+k:T,:)' * y_demeaned(1:T-k,:));
end

% Asymptotic bands
critval = norminv(1-alph);
ul = repmat(critval/sqrt(T),pmax,1);
ll = -1*ul;

% Barplots
figure('name','Autocorrelation');
bar(r_k);
hold on;
plot(1:pmax,ul,'color','black','linestyle','--');
plot(1:pmax,ll,'color','black','linestyle','--');
hold off;

% The following is just for pretty plots
acfbarplot = gca; % Get current axes handle
acfbarplot.Title.String = 'Sample autocorrelation coefficients';
acfbarplot.XAxis.Label.String = 'lags';
acfbarplot.XAxis.TickValues = 1:pmax;
acfbarplot.YAxis.Label.String = 'acf value';
acfbarplot.YAxis.Limits = [-1 1];
acfbarplot.XAxis.Limits = [0 pmax];

end % function end