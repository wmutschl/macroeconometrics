% =========================================================================
% Declare model equations for:
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
model;
% auxiliary expressions, these will be substituted by Dynare preprocessor
@#if CRRA
    #UC  = gam*C^(-etaC);
    #UCp = gam*C(+1)^(-etaC);
    #UL  = -pssi*(1-L)^(-etaL);
@# else 
    #UC  = gam*C^(-1);
    #UCp = gam*C(+1)^(-1);
    #UL  = -pssi*(1-L)^(-1);
@#endif

[name='Euler equation']
UC=betta*UCp*(1-delt+R);

[name='labor supply']
W=-UL/UC;

[name='capital law of motion'] 
K=(1-delt)*K(-1)+I;

[name='market clearing/resource constraint']
Y=I+C;

[name='production function']
Y=A*K(-1)^alph*L^(1-alph);

[name='labor demand']
W=(1-alph)*Y/L;

[name='capital demand']
R=alph*Y/K(-1);

[name='exogenous TFP process']
log(A)=rhoA*log(A(-1))+eps_A;

end;