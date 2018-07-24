% Set up model
A = [0.2, 0.3; 
    -0.6, 1.1];
SIGu = [0.9, 0.2;
        0.2, 0.5];
nA = size(A,1);
% Theoretical Moments VAR(1)
mu = zeros(2,1);
Gamma0 = dlyapdoubling(A,SIGu);

% Simulate data
R = 10000; 
T = 1000;
% Intialize storage for Monte-Carlo
meanVAR = nan(nA,1,R);
covVAR = nan(nA,nA,R);
waitb = waitbar(0,'Number of MC runs');
for r = 1:R    
    % Intialize storage for observations in Monte Carlo run r
    % First observation are zeros, ie equal to expectation
    y = zeros(2,T);
    u = transpose(mvnrnd(mu,SIGu,T)); % draw T times from multivariate normal
    for t=2:T
        y(:,t) = A*y(:,t-1) + u(:,t);
    end
    % Get sample statistics for Monte Carlo sample r
    % Either use build-in functions or formulas given in lecture
    meanVAR(:,:,r) = mean(y,2);
    covVAR(:,:,r) = cov(y');
    waitbar(r/R);
end
close(waitb);
%Compare theoretical with simulated moments for VAR(1)
[mean(meanVAR,3) mu]
[mean(covVAR,3) Gamma0]