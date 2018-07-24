% =======================================================================
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
% A and B are the coefficient matrices of the difference equation
% A*z(t+1) = B*z(t)
% where z(t) is arranged so that the state variables come first, and
% nx is the number of state variables.%
% Outputs: the decision rule hx and the law of motion gx. If we write
% z(t) = [x(t);y(t)] where x(t) contains precisely the state variables, then
% y(t) = gx*x(t) and
% x(t+1) = hx*x(t) +sig*eta*epsilon(t+1)
% =======================================================================
% Willi Mutschler, January 2018
% willi@mutschler.eu
% =======================================================================
%% Symbolic definitions
%Define the structural parameters
syms alph betta rhoA sigA
%Define the state variables in this period (_cu) and the next period(_cup)
syms K_cu K_cup A_cu A_cup
%Define the control variables in this period (_cu) and the next period (_cup)
syms C_cu C_cup

% Define the vector of states, x and xp
x  = [K_cu  A_cu ];
xp = [K_cup A_cup];
% Define the vector of controls, y and yp
y =  [C_cu];
yp = [C_cup];

nx  = size(x,2);
ny  = size(y,2);
ne = 1;
n   = nx+ny;
% The eta matrix: nx*ne
eta      = sym(zeros(nx,ne));
eta(2,1) = sigA;

%% Model equations 
% Eq 1: Euler
eq1 = -C_cu^(-1) + alph*betta*C_cup^(-1)*A_cup*K_cup^(alph-1);
% Eq 2: Capital law of motion
eq2 = -K_cup + A_cu*K_cu^alph - C_cu;
% Eq 3; Exogenous TFP process
eq3 = -log(A_cup) + rhoA*log(A_cu);

% Create function f
f = [eq1;eq2;eq3];

%% steady state
Abar = sym(1);
Kbar = (alph*betta*Abar)^(1/(1-alph));
Cbar = Abar*Kbar^alph-Kbar;
xbar = [Kbar Abar];
ybar = Cbar;

%% The first order derivatives
f1 = jacobian(f,xp);
f2 = jacobian(f,yp);
f3  = jacobian(f,x);
f4  = jacobian(f,y);

%% Set numerical values and evaluate steady state
alph  = 0.35;
betta = 0.99;
rhoA  = 0.9;
sigA  = 0.6;
K_cu = eval(subs(Kbar)); K_cup = K_cu;
A_cu = eval(subs(Abar)); A_cup = A_cu;
C_cu = eval(subs(Cbar)); C_cup = C_cu;
nf   = eval(subs(f)); 
nf1 = eval(subs(f1)); 
nf2 = eval(subs(f2)); 
nf3  = eval(subs(f3)); 
nf4  = eval(subs(f4)); 
neta = eval(subs(eta)); 



%% Approximation to first order
A = [-nf1 -nf2];
B = [nf3 nf4];
realsmall=1e-7;
try
    [S,T,Q,Z] = qz(A,B); % upper triangular factorization of the matrix pencil
catch
    error('Solab:noexist','Error using qz');
end
disp(Q*A*Z-S)
disp(Q*B*Z-T)
% Rule out that both s_ii and t_ii are zero
zxz = sum((abs(diag(S))<realsmall) & (abs(diag(T))<realsmall));
if ~(~zxz)
	error('Solab:noexist','Coincident zeros');
end
%   If S is triangular, the diagonal elements of S and T,
%       sii = diag(S), tii = diag(T),
%   are the generalized eigenvalues that satisfy
%       A*V*diag(tii) = B*V*diag(sii)
%       diag(tii)*W'*A = diag(sii)*W'*B
%   where matrices V's and W' columns are the generalized eigenvectors
%   The eigenvalues produced by
%       lambda = eig(A,B)
%   are the ratios of the sii and tii.
%       lambda = sii./tii
eig(A,B)
disp(diag(S)./diag(T));
% reordering such that stable (smaller than one) generalized Eigenvalues 
% of B w.r.t. A is in the upper left of S and T, abs(tii)>abs(sii)
[S,T,Q,Z] = ordqz(S,T,Q,Z,'udo');
disp(abs(diag(S))./abs(diag(T)));
% Blanchard-Khan order condition:
if abs(T(nx,nx))>abs(S(nx,nx))
	error('Solab:noexist','No equilibrium exists.');
elseif abs(T(nx+1,nx+1))<abs(S(nx+1,nx+1))
	error('Solab:indeterminacy','Indeterminacy.');    
end
z21 = Z(nx+1:end,1:nx);
z11 = Z(1:nx,1:nx);
% Blanchard-Khan rank condition:
if rank(z11)<nx
    error('Solab:noexist','Invertibility condition violated')
end
z11i = z11\eye(nx);
s11 = S(1:nx,1:nx);
t11 = T(1:nx,1:nx);

dyn = s11\t11;
% real function takes away very small imaginary parts of the solution
gx = real(z21*z11i); 
hx = real(z11*dyn*z11i);

%% Equivalence to Dynare's policy and transition function
% note that in Dynare K_cu = K(-1) and A_cu = A(-1)
PolFct = [eval(Cbar)      eval(Kbar)        eval(Abar);
          gx(1)           hx(1,1)           hx(2,1);
          gx(2)*hx(2,2)   hx(1,2)*hx(2,2)   hx(2,2);
          gx(2)*neta(2,1) hx(1,2)*neta(2,1) neta(2,1)];
rowNames = {'Constant','K(-1)','A(-1)','eps_A'};
colNames = {'C','K','A'};
PolFctDisp = array2table(PolFct,'RowNames',rowNames,'VariableNames',colNames);
disp(PolFctDisp)