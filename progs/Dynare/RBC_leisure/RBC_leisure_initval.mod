% =========================================================================
% Define initial values for steady state computation for:
% -------------------------------------------------------------------------
% Baseline RBC model with 
% - Consumption and leisure in either CRRA or log utility function, 
%   dependent on Dynare macro CRRA.
% - Cobb Douglas production function with labor & capital in a perfect
%   competition setting
% - Total Factor Productivity shock
% -------------------------------------------------------------------------
% Based on codes by Martin Andreasen and Johannes Pfeifer
% Willi Mutschler, February 2018
% willi@mutschler.eu
% =========================================================================
initval;
    A =  1;
    Y =  1;
    C =  0.8;
    K =  10;
    L =  1/3;
    R =  0.03;
    W =  2.2;
    I =  0.3;
end;