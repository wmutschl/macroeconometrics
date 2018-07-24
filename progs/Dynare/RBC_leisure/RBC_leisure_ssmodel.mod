% =========================================================================
% Steady_state_model block for:
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in log utility function (not CRRA)
% - Cobb Douglas production function with labor & capital in a perfect
%   competition setting
% - Total Factor Productivity shock
% -------------------------------------------------------------------------
% Based on codes by Martin Andreasen and Johannes Pfeifer
% Willi Mutschler, February 2018
% willi@mutschler.eu
% =========================================================================
steady_state_model;
    A = 1;                             % technology innovations
    R = 1/betta+delt-1;                % rental rate of capital
    K_L = ((alph*A)/R)^(1/(1-alph));   % capital divided by labor
    W = (1-alph)*A*K_L^alph;           % wage level
    I_L = delt*K_L;                    % investment divided by labor
    Y_L = A*K_L^alph;                  % output divided by labor
    C_L = Y_L-I_L;                     % consumption divided by labor
    % closed-form expression for labor when using log utility
    L = gam/pssi*C_L^(-1)*W/(1+gam/pssi*C_L^(-1)*W);
    C = C_L*L; % consumption level
    I = I_L*L; % investment level
    K = K_L*L; % capital level
    Y = Y_L*L; % output level
end;
