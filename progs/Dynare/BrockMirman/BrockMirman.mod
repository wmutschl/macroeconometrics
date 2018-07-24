% =========================================================================
% Stochastic growth model of Brock and Mirman (1972) with technology shock
% Computes simulated data and impulse response functions based on the true 
% decision functions and compares these with the corresponding objects 
% based on Dynare's approximated model solutions
% =========================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =========================================================================

% =========================================================================
% Set options for Dynare preprocessor
% =========================================================================
% comparison: 1 for comparing simulated data and IRFs based on true
%             solution vs. Dynare's approximated one
@#define comparision = 0

% =========================================================================
% Declare variables, shocks and parameters
% =========================================================================
var 
    C        ${C}$ (long_name='consumption')
    K        ${K}$ (long_name='capital')
    A        ${Z}$ (long_name='total factor productivity')
;

varexo
    eps_A    ${\varepsilon_A}$ (long_name='TFP shock')
;
    
parameters
    alph    ${\alpha}$   (long_name='capital share')
    betta   ${\beta}$    (long_name='discount factor')
    rhoA    ${\rho_A}$   (long_name='persistence TFP')
    sigA    ${\sigma_A}$ (long_name='standard deviation TFP shock')
;

% =========================================================================
% Calibrate parameter values based on OECD data
% =========================================================================
alph  = 0.35;
betta = 0.99;
rhoA  = 0.9;
sigA  = 0.6;

% =========================================================================
% Model equations
% =========================================================================
model;
[name='Euler equation']
C^(-1)=alph*betta*C(+1)^(-1)*A(+1)*K^(alph-1);
[name='capital law of motion'] 
K=A*K(-1)^alph-C;
[name='exogenous TFP process']
log(A)=rhoA*log(A(-1))+sigA*eps_A;
end;

% =========================================================================
% Define shock covariance matrix
% =========================================================================
shocks;
    var eps_A = 1;
end;

% =========================================================================
% Define steady state computation
% =========================================================================
% Analytical steady state using steady_state_model
steady_state_model;
    A = 1;                           % technology level
    K = (alph*betta*A)^(1/(1-alph)); % capital level
    C = A*K^alph-K;                  % consumption level
end;

% =========================================================================
% Computations
% =========================================================================
steady; % compute steady state given the starting values
resid;  % check the starting values for the steady state
check;  % check Blanchard & Khan rank condition
@# if comparision == 0
    stoch_simul(order=1);
@# else
    stoch_simul(order=1,irf=40,periods=200);
    % order=1:     first-order approximation of solution
    % irf=40 :     compute IRFs for 40 periods for one-standard deviation shock
    %              from approximated solution
    % periods=200: draw 200 shocks and simulate 200 observations 
    %              from approximated solution

    %% Compare Trajectories
    % This is the function used by DYNARE to simulate 200 observations from the
    % approximated solution, note that these are already saved in C K A
    epsA = oo_.exo_simul;
    ybar = oo_.steady_state;
    yapprox=simult_(ybar,oo_.dr,epsA,options_.order);
    isequal(yapprox(:,2:end),[C K A]')

    % Simulate the model using the exact policy functions given the same draws 
    % of exogenous shocks as was used by DYNARE
    ytrue = zeros(size(ybar,1),options_.periods+1);
    ytrue(:,1) = ybar;
    for i = 2:(options_.periods+1)
        ytrue(3,i) = ytrue(3,i-1)^rhoA * exp(epsA(i-1));           % technology
        ytrue(2,i) = alph*betta*ytrue(3,i)*ytrue(2,i-1)^alph;      % capital
        ytrue(1,i) = (1-alph*betta)*ytrue(3,i)*ytrue(2,i-1)^alph;  % consumption
    end

    % Plot trajectories of exact and approximated solution
    figure('Name','Trajectories of Simulated Data')
    subplot(3,1,1);
        plot((ytrue(1,2:end)));
        hold on
        title('Consumption');
        plot(C-ybar(1),'--')
        hlineC = refline([0 ybar(3)]);
        hlineC.Color = 'black';
        hlineC.LineStyle = ':';
        xlim([1 options_.periods]);
        hold off
        legend('Exact Solution','Approximated Solution','Steady-State')
    subplot(3,1,2);
        plot(log(ytrue(2,2:end)));
        hold on
        title('Capital');
        plot(K-ybar(2),'--')
        hlineK = refline([0 ybar(3)]);
        hlineK.Color = 'black';
        hlineK.LineStyle = ':';
        xlim([1 options_.periods]);
        hold off
        legend('Exact Solution','Approximated Solution','Steady-State')
    subplot(3,1,3);
        plot(log(ytrue(3,2:end)));
        hold on
        title('Technology');
        plot(A-ybar(3),'--')
        hlineA = refline([0 ybar(3)]);
        hlineA.Color = 'black';
        hlineA.LineStyle = ':';
        xlim([1 options_.periods]);
        hold off
        legend('Exact Solution','Approximated Solution','Steady-State')    


    %% Compare Impulse Response Functions
    % This is the function used by DYNARE to compute IRFs from the approximated
    % solution, note that these are already saved in C_eps_A, K_eps_A A_epsA
    onestddev = chol(M_.Sigma_e+1e-14*eye(M_.exo_nbr));
    irfapprox=irf(oo_.dr,onestddev, options_.irf, options_.drop, options_.replic, options_.order);
    isequal(irfapprox,[C_eps_A K_eps_A A_eps_A]')

    % Compute IRFs using the exact policy functions
    irftrue = zeros(size(ybar,1),options_.periods+1);
    epsA2 = zeros(options_.irf,M_.exo_nbr);
    epsA2(1,:) = onestddev;
    irftrue = zeros(size(ybar,1),options_.irf+1);
    irftrue(:,1) = ybar;
    for i = 2:(options_.irf+1)
        irftrue(3,i) = irftrue(3,i-1)^rhoA * exp(epsA2(i-1));    % technology
        irftrue(2,i) = alph*betta*irftrue(3,i)*irftrue(2,i-1)^alph;      % capital
        irftrue(1,i) = (1-alph*betta)*irftrue(3,i)*irftrue(2,i-1)^alph;  % consumption
    end    

    % Plot IRFs of exact and approximated solution
    figure('Name','Impulse Response Functions');
    subplot(3,1,1);
        plot(log(irftrue(1,2:end))-log(ybar(1)));
        hold on
        title('Consumption');
        plot(C_eps_A,'--')
        hlineC = refline([0 0]);
        hlineC.Color = 'black';
        hlineC.LineStyle = ':';
        xlim([1 options_.irf]);
        hold off
        legend('Exact IRF','Approximated IRF','Steady-State')    
    subplot(3,1,2);
        plot(log(irftrue(2,2:end))-log(ybar(2)));
        hold on
        title('Capital');
        plot(K_eps_A,'--')
        hlineK = refline([0 0]);
        hlineK.Color = 'black';
        hlineK.LineStyle = ':';
        xlim([1 options_.irf]);
        hold off
        legend('Exact IRF','Approximated IRF','Steady-State')    
    subplot(3,1,3);
        plot(log(irftrue(3,2:end))-log(ybar(3)));
        hold on
        title('Technology');
        plot(A_eps_A,'--')
        hlineA = refline([0 0]);
        hlineA.Color = 'black';
        hlineA.LineStyle = ':';
        xlim([1 options_.irf]);
        hold off
        legend('Exact IRF','Approximated IRF','Steady-State')
@# endif