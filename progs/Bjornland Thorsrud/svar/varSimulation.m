function [alfab,betab,sigmab]=varSimulation(alfaC,betaC,emat,nlag,draws)
% PURPOSE: Simulate a VAR DGP by bootstrapping, e.g., by
% drawing from the residuals. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   betaC = matrix. (n+(n*(nlag-1))) x (m-1)). VAR(p) respresented as
%   VAR(1). beta parameters in VAR regression (not including const.)
%
%   alfaC = vector. (n+(n*(nlag-1))) x 1). VAR(p) respresented as
%   VAR(1). Constant in VAR regression
%
%   emat = residual matrix. (t x n), where t is the number of observations
%   and n is the number of variables in the VAR
%   
%   nlag = number of lags used in VAR
%   
% Output:
%   
%   betab = matrix. This will be the bootstrapped beta parameter values. 
%   (nvar x nvar*nlag x draws)
%
%   alfab = matrix. This will be the bootstrapped alfa parameter values. 
%   (nvar x 1 x draws)
%
%   sigmab = covariance matrix from each draw, (nvar x nvar x draws)
%       
%   IMPORTANT: The variables needs to be ordered the following way in the 
%   betaC coeff matrix:
%   [y1(t-1) y2(t-1) yn(t-1)...y1(t-nlag) y2(t-nlag) yn(t-nlag)]
%   The alfaC and betaC vector and matrix must also represent VAR(p) as
%   VAR(1), see function varGetCompForm.m
%
% Usage:
% [alfab,betab,sigmab]=varSimulation(alfaC,betaC,emat,nvar,nlag,T,draws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% 2010
% leifath@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,nvar]=size(emat);

% Some initial stuff...
options.draws=draws;
options.burninn=10;
options.tPlus=0;
options.stability=true;
if options.tPlus>options.burninn-30
    options.burninn=2*options.tPlus;
end;
if isempty(alfaC)
    options.constant=false;
    alfaC=zeros(size(betaC,1),1);
else
    options.constant=true;
end;

% Empty output
alfab=nan(nvar,1,options.draws);
betab=nan(nvar,nvar*nlag,options.draws);
sigmab=nan(nvar,nvar,options.draws);

zeroAddOn=zeros((nlag-1)*nvar,1);
ematm=emat-repmat(mean(emat,1),[T 1]);        

for d=1:options.draws    
    
    stabilityCondition=1;
    
    while stabilityCondition==1        
        
        
        drawsIdx=ceil(rand((T+options.burninn),1)*T);
        randResid=ematm(drawsIdx,:)';                                                

        %empty y mat for iterated vector
        y=nan(nvar,options.burninn+T);
        yy=zeros(nvar*nlag,1);
        
        %simulate system                
        for s=1:options.burninn+T
            ee=[randResid(:,s);zeroAddOn];            
            yy=alfaC+betaC*yy+ee;            
            y(:,s)=yy(1:nvar,1);
        end;
        
        % estimate new parameters         
        xl=latMlag(y',nlag);
        if options.constant
            xl=[ones(T,1) xl(end-T+1:end,:)];
        else
            xl=xl(end-T+1:end,:);
        end;
        yl=y(:,end-T+1:end)';
        
        aabb=(xl\yl)';
        resid=yl-(aabb*xl')';           
        sigma=(resid'*resid)/T;
                           
        % get comp. form to check stability
        if options.stability
            if options.constant
                betaCbb=varGetCompForm(aabb(:,2:end),[],nlag,nvar);
            else
                betaCbb=varGetCompForm(aabb,[],nlag,nvar);
            end;
            if any(abs(real(eig(betaCbb)))>=1)
                stabilityCondition=1;
            else
                % To brake out of while statement if stability
                % ensured
                stabilityCondition=0;
            end;
        else            
            % To brake out of while statement if stability not important
            stabilityCondition=0;
        end;    
    end;
    
    if options.constant
        alfab(:,1,d)=aabb(:,1);
        betab(:,:,d)=aabb(:,2:end);
    else
        betab(:,:,d)=aabb;        
    end;    
    sigmab(:,:,d)=sigma;             
    
end;

