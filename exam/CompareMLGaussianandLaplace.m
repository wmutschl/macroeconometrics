clear, clc, close all

% Simulate data for AR(1) with c=1, phi=0.8 and Laplace distribution
c = 1;
phi = 0.8;
sigm = sqrt(2);
T = 200;
y = nan(T,1);
y(1) = c/(1-phi);
for t=2:T
    y(t) = c + phi*y(t-1)+sigm*randl(1,1);
end

size(y)
save('LaPlace.txt', 'y','-ascii','-tabs');
clear 
load LaPlace.txt
%y = xlsread('../data/data.xlsx','AR4');
% True values are c=1; phi=[0.51,-0.1, 0.06, -0.22] and sigma=0.5; 
p = 1;
const = 1;
ML0  = ARpML(LaPlace,p,const,0);
ML1  = ARpML(LaPlace,p,const,1);
[ML0.thetatilde ML1.thetatilde]