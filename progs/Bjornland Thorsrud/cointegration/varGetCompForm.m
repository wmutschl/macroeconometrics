function [BETAC,ALFAC]=varGetCompForm(beta,alfa,nlag,nvar)
% PURPOSE: Represent VAR(p) as VAR(1), eg comp form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   alfa = vector. (n x 1), where n is the number of variables. Constants
%   form regression
%
%   beta = matrix. (n x (m-1)), where n is the number of variables, and m-1
%   is the number of variables*number of lags. 
%
%   nlag = number of lags in VAR 
%
%   nvar = number of variables in VAR
%
% Output:
%   betaC = matrix. (n+(n*(nlag-1))) x (m-1)). 
%
%   alfaC = vector. (n+(n*(nlag-1))) x 1). 
%
% Usage:
%   [BETAC,ALFAC]=varGetCompForm(beta,alfa,nlag,nvar) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leif Anders Thorsrud
% 2010
% leifath@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE:
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,m1,d]=size(beta);
BETAC=nan(m1,m1,d);
if ~isempty(alfa)
    ALFAC=nan(m1,1,d);
else
    ALFAC=[];
end

for i=1:d
    if nlag>1
        BETAC(:,:,i)=[beta(:,:,i);
            eye((nvar)*(nlag-1)) zeros((nvar)*(nlag-1),nvar)];    
        if ~isempty(alfa)
            ALFAC(:,:,i)=[alfa(:,i);    
                zeros((nvar)*(nlag-1),1)];        
        end;
    else
        BETAC(:,:,i)=beta(:,:,i);
        if ~isempty(alfa)
            ALFAC(:,:,i)=alfa(:,i);    
        end;            
    end;                
end;
