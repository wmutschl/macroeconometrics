clear, clc, close all

% Simulate data for AR(1) with c=1, phi=0.8 and Laplace distribution
c = 1;
phi = 0.8;
sigm = sqrt(2);
T = 1000;
y = nan(T,1);
y(1) = c/(1-phi);
for t=2:T
    y(t) = c + phi*y(t-1)+sigm*randl(1,1);
end

y = y(801:900);
size(y)