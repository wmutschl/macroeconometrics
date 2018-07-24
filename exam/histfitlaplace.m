%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Histogram with superimposed fitted        %
% Laplace distribution with MATLAB Implementation %
%                                                 %
% Author: M.Sc. Eng. Hristo Zhivomirov   04/29/15 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function histfitlaplace(x)

% function: histfitlaplace(x)
% x - data sequence

% distribution parameters estimation
mu = mean(x);                                                   % find the distribution expected value
sigma = std(x);                                                 % find the distribution std
xprim = mu + linspace(-3*sigma, 3*sigma, 1000);                 % xprim vector generation [-3*sigma+mu 3*sigma+mu]
fx = 1/sqrt(2*(sigma^2))*exp(-sqrt(2/(sigma^2))*abs(xprim-mu)); % calculate the Laplace PDF

% normalize the pdf using the area of the histogram
bins = sqrt(length(x)); % number of the histogram bins
n = length(x);          % length of the data
r = max(x) - min(x);    % range of the data
binwidth = r/bins;      % width of each bin
area = n*binwidth;      % area of the histogram
fx = area*fx;           % normalize the pdf

% plot the histogram
hist(x, bins);
grid on
hold on

% plot the Laplace PDF 
plot(xprim, fx, '-r', 'LineWidth', 2)

% set the x limits
dxl = mu - min(x);
dxr = max(x) - mu;
dx = max([dxl dxr]);
xlim([mu-1.1*dx mu+1.1*dx])

end