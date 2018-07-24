function [a,P,LL,ap,Pp]=KalmanFilter(y,Z,d,H,T,c,R,Q,a0,P0)
% PURPOSE:	Apply the Kalman Filter to a time-varying state space
% formulation, specified as in Harvey (1989). Faster version of 
% KalmanFilterOrig, i.e no time variation allowed, no error checking, and
% no analytic derivation option for derivatives. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage: 
%
% [a,P,LL,ap,Pp]=KalmanFilter(y,Z,d,H,T,c,R,Q,a0,P0)
%


% 1. Initialisation
%[num_observables,num_observations]=size(y);
num_observations=size(y,2);
num_states=size(a0,1);
% Pre-allocate space for state and covariance of state
[a,ap]=deal(NaN(num_states,1,num_observations));           % state / conditional 
[P,Pp]=deal(NaN(num_states,num_states,num_observations));  % covariance of state / conditional
LL=zeros(1,num_observations);

% 2. The first observation    
ap(:,1,1)=T*a0+c;
Pp(:,:,1)=T*P0*T'+R*Q*R';

mask=~isnan(y(:,1));
W=diag(mask);
Wt=W(mask,:);
    
Ft=Wt*(Z*Pp(:,:,1)*Z'+H)*Wt';
Ft_inverse=inv(Ft);
y(~mask,1)=0;
vt = Wt*(y(:,1)-Z*ap(:,1,1)-d);
    
a(:,1,1)=ap(:,1,1)+Pp(:,:,1)*Z'*Wt'*Ft_inverse*vt;
P(:,:,1)=Pp(:,:,1)-Pp(:,:,1)*Z'*Wt'*Ft_inverse*Wt*Z*Pp(:,:,1);
                         
% 3. Subsequent observations
for t=2:1:num_observations;
    
    ap(:,1,t)=T*a(:,t-1)+c;
    Pp(:,:,t)=T*P(:,:,t-1)*T'+R*Q*R';     
    
    mask=~isnan(y(:,t));
    W=diag(mask);
    Wt=W(mask,:);
    
    if ~isempty(Wt); % we have some current data 
        Ft = Wt*(Z*Pp(:,:,t)*Z'+H)*Wt';    % Variance of observables
        
        Ft_inverse=inv(Ft);
        y(~mask,t)=0; %replaces nan's with zeros
        vt = Wt*(y(:,t)-Z*ap(:,1,t)-d);
      
        a(:,1,t)=ap(:,1,t)+Pp(:,:,t)*Z'*Wt'*Ft_inverse*vt;               
        P(:,:,t)=(Pp(:,:,t)-Pp(:,:,t)*Z'*Wt'*Ft_inverse*Wt*Z*Pp(:,:,t));    
        
        % add likelihood here if needed
        LL(:,t)=-0.5*((sum(mask))'*log(2*pi)+log(det(Ft))+vt'*Ft_inverse*vt);                    
        
    else % we have no data for the current period (so our expectations are unchanged)      
        a(:,:,t)=ap(:,1,t);
        P(:,:,t)=Pp(:,:,t);
    end;
end;




