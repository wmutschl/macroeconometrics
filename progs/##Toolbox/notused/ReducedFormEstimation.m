function [Ahat,Uhat,SigmahatU,SigmaAhat] = ReducedFormEstimation(Z,Y,EstimationMethod,Model)

T = size(Y,2);

% Reduced-form estimation
if strcmp(Model,'LEVELS')
    switch EstimationMethod
        case 1 %Unrestricted Least Squares 
            % OLS
            %Ahat = Y*Z'*inv(Z*transpose(Z)); %using svd
            %Uhat = Y-Ahat*Z;        
            [AhatT,UhatT] = ols(Y',Z'); %using svd
            Ahat = transpose(AhatT); Uhat = transpose(UhatT);        
            SigmahatU = Uhat*transpose(Uhat)/(T-size(Z,1));
            SigmaAhat = kron((Z*transpose(Z))\eye(size(Z,1)),SigmahatU);
        case 4 % Gaussian ML        
            % Gaussian ML is equivalent to OLS, except SigmahatU
            [AhatT,UhatT] = ols(Y',Z'); %using svd
            Ahat = transpose(AhatT); Uhat = transpose(UhatT);        
            SigmahatU = Uhat*transpose(Uhat)/(T);
            SigmaAhat = kron((Z*transpose(Z))\eye(size(Z,1)),SigmahatU);
    end
elseif strcmp(Model,'VECM')
end