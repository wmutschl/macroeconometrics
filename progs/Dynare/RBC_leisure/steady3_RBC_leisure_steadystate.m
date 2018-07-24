function [ys,check] = RBC_leisure_steadystate(ys_, exo_)
% function [ys,check] = RBC_leisure_steadystate(ys,exo)
% Computes the steady state for the RBC_LEISURE_CRRA.mod and uses a numerical
% solver to do so
% Inputs: 
%   - ys        [vector] vector of initial values for the steady state of
%                   the endogenous variables
%   - exo       [vector] vector of values for the exogenous variables
%
% Output: 
%   - ys        [vector] vector of steady state values for the the endogenous variables
%   - check     [scalar] set to 0 if steady state computation worked and to
%                    1 if not (allows to impose restriction on parameters)
% Based on codes by Martin Andreasen
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================


global M_    % Get Dynare global model structure

%% read out parameters to access them with their name
NumberOfParameters = M_.param_nbr;
for ii = 1:NumberOfParameters
  paramname = deblank(M_.param_names(ii,:));
  eval([ paramname ' = M_.params(' int2str(ii) ');']);
end

%% Intializations
check = 0;          % initialize indicator: 0 for no errros, else 1
options=optimset(); % set options for numerical solver

%% own steady state computation
dy = 0;                          % output growth
A = 1;                           % technology innovations
R = 1/betta+delt-1;              % rental rate of capital
K_L = ((alph*A)/R)^(1/(1-alph)); % capital divided by labor
if K_L <= 0
    check = 1; % set failure indicator
    return;    % return without updating steady states
end
W = (1-alph)*A*K_L^alph;         % wage level
I_L = delt*K_L;                   % investment divided by labor
Y_L = A*K_L^alph;                 % output divided by labor
C_L = Y_L-I_L;                     % consumption divided by labor
if C_L <= 0
    check = 1; % set failure indicator
    return;    % return without updating steady states
end

% The labor level
if etaC == 1 && etaL == 1
    % Closed-form solution for labor
    L = gam/pssi*C_L^(-1)*W/(1+gam/pssi*C_L^(-1)*W);
else
    % No closed-form solution use a fixed-point algorithm
    L0 = 1/3;
    [L,~,exitflag] = fsolve(@(L) pssi*(1-L)^(-etaL)*L^etaC - gam*C_L^(-etaC)*W, L0,options);    
    if exitflag <= 0
        check = 1; % set failure indicator
        return     % return without updating steady states
    end    
end

% Value of remaining variables
C = C_L*L; % consumption level
I = I_L*L; % investment level
K = K_L*L; % capital level
Y = Y_L*L; % output level

%% write out variables independent of name
NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = deblank(M_.endo_names(ii,:));
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
ys = ys';