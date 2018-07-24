function [smoothState, smoothVarCov]=KalmanSmoother(y,a,P,ap,Pp,T)
% PURPOSE:	Apply the Kalman Smoother to a state space formulation as in Harvey (1989)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [~,num_observations]=size(y);   
    
 as(:,:,num_observations)=a(:,:,num_observations);
 Ps(:,:,num_observations)=P(:,:,num_observations);
 
 for t=num_observations-1:-1:1;
     Pt=P(:,:,t);
     Ppt_lead=Pp(:,:,t+1);
     
     Pstar=Pt*T'*pinv(Ppt_lead);    
     Ps(:,:,t)=Pt+Pstar*(Ps(:,:,t+1)-Ppt_lead)*Pstar';
    
     as(:,:,t)=a(:,:,t)+Pstar*(as(:,:,t+1)-ap(:,:,t+1));   
 end;
 smoothState=as;
 smoothVarCov=Ps;
 
 
 