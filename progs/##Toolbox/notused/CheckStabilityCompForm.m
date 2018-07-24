function Stability = CheckStability(maxEig)
% Stability of var
if maxEig>=1
    error('CheckStabilityCompForm:Stability','The VAR is not stationary!')
    Stability = 0;
else
    disp('The VAR is stationary!');
    Stability = 1;
end;
