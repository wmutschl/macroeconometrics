% =======================================================================
% Illustration of the weak law of large numbers for several distributions 
% and AR(1) process based on these distributions. Distributions considered:
% normal,uniform,geometric,student's (finite and infinite variance),gamma.
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================
clearvars; clc;close all;
% Initializations
T = 10000;                                  % maximum horizon of periods
z = nan(1,T); sig_z = 2; mu_z = 10;         % normal distribution
u = nan(1,T); a = 2; b = 4; mu_u = (a+b)/2; % uniform distribution
ge = nan(1,T); p = 0.2; mu_ge = (1-p)/p;    % geometric distribution
ga = nan(1,T); k=2; thet=2; mu_ga=k*thet;   % gamma distribution
st1 = nan(1,T); nu1 = 3; mu_st1 = 0;        % student t with finite variance
st2 = nan(1,T); nu2 = 2; mu_st2 = 0;        % student t with infinite variance
Y = nan(6,T); phi=0.8; mu_y = 0;            % stable AR(1) process

% Draw random variables
ZZ = randn(1,T);                            % normal distribution
UU = rand(1,T);                             % uniform distribution
GeGe = geornd(p,1,T);                       % geometric distribution
GaGa = gamrnd(k,thet,1,T);                  % gamma distribution
ST1 = trnd(nu1,1,T);                        % student t with finite variance
ST2 = trnd(nu2,1,T);                        % student t with infinite variance

wait = waitbar(0,'Please wait...'); % open waitbar
for t = 1:T
    Zt = mu_z + sig_z.*ZZ(1:t);             % normal distribution
    Ut = (b-a).*UU(1:t) + a;                % uniform distribution
    Get = GeGe(1:t);                        % geometric distribution
    Gat = GaGa(1:t);                        % gamma distribution
    ST1t = ST1(1:t);                        % student t with finite variance
    ST2t = ST2(1:t);                        % student t with infinite variance
    Yt(:,1) = zeros(6,1);                   % initialize AR(1) processes
    if t>1
        for tt=2:t
            % Note that we demean the errors
            Yt(1,t) = phi*Yt(1,t-1) + (Zt(t)-mu_z);     % normal
            Yt(2,t) = phi*Yt(2,t-1) + (Ut(t)-mu_u);     % uniform
            Yt(3,t) = phi*Yt(3,t-1) + (Get(t)-mu_ge);   % geometric
            Yt(4,t) = phi*Yt(4,t-1) + (Gat(t)-mu_ga);   % gamma 
            Yt(5,t) = phi*Yt(5,t-1) + (ST1t(t)-mu_st1); % finite student t
            Yt(6,t) = phi*Yt(6,t-1) + (ST2t(t)-mu_st2); % infinite student t
        end
    end
    % Compute and store averages
    z(t) = mean(Zt);
    u(t) = mean(Ut);
    ge(t) = mean(Get);
    ga(t) = mean(Gat);
    st1(t) = mean(ST1t);
    st2(t) = mean(ST2t);
    Y(:,t) = mean(Yt,2);
    waitbar(t/T); % update waitbar
end
close(wait); % close waitbar

% Plots for distributions
figure('name','Law of Large Numbers for different distributions');
subplot(2,3,1); 
    plot(z); 
    line(0:T,repmat(mu_z,1,T+1),'linestyle','--','color','black'); 
    title('Normal');
subplot(2,3,2); 
    plot(u); 
    line(0:T,repmat(mu_u,1,T+1),'linestyle','--','color','black');
    title('Uniform');
subplot(2,3,3); 
    plot(ge); 
    line(0:T,repmat(mu_ge,1,T+1),'linestyle','--','color','black'); 
    title('Geometric');
subplot(2,3,4); 
    plot(ga); 
    line(0:T,repmat(mu_ga,1,T+1),'linestyle','--','color','black'); 
    title('Gamma');
subplot(2,3,5); 
    plot(st1); 
    line(0:T,repmat(mu_st1,1,T+1),'linestyle','--','color','black'); 
    title('Student''s t finite variance');
subplot(2,3,6); 
    plot(st2); 
    line(0:T,repmat(mu_st2,1,T+1),'linestyle','--','color','black'); 
    title('Student''s t infinite variance');

% Plots for AR(1) processes
figure('name','Law of Large Numbers for stable AR(1)');
subplot(2,3,1); 
    plot(Y(1,:)); 
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black'); 
    title('Normal errors');
subplot(2,3,2); 
    plot(Y(2,:)); 
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black');
    title('Uniform errors');
subplot(2,3,3); 
    plot(Y(3,:)); 
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black'); 
    title('Geometric errors');
subplot(2,3,4); 
    plot(Y(4,:));  
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black'); 
    title('Gamma errors');
subplot(2,3,5); 
    plot(Y(5,:));  
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black'); 
    title('Student''s t finite variance errors');
subplot(2,3,6); 
    plot(Y(6,:));
    line(0:T,repmat(mu_y,1,T+1),'linestyle','--','color','black'); 
    title('Student''s t infinite variance errors');        