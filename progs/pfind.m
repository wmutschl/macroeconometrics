% PFIND.M
% Lutz Kilian
% University of Michigan
% November 2008

function [sichat,hqchat,aichat]=pfind(y,pmax);

[t,K]=size(y);

% Construct regressor matrix and dependent variable
XMAX=ones(1,t-pmax);
for i=1:pmax
	XMAX=[XMAX; y(pmax+1-i:t-i,:)'];
end
Y=y(pmax+1:t,:)';   

% Evaluate criterion for p=0,...,pmax
for jj=0:pmax

	X=XMAX(1:jj*K+1,:);

	B=Y*X'*inv(X*X');        
    SIGMA=(Y-B*X)*(Y-B*X)'/(t-pmax);
	np=length(B(:));      % Number of freely estimated parameters
    
   % Lutkepohl suggests np=length(vec(B))-K which is used in table 4.5.
   % This does not affect the ranking, but the value of the criterion
   % function. See p. 147.
   
    aiccrit(jj+1,1)=log(det(SIGMA))+np*2/(t-pmax);    	   	         % AIC value
    hqccrit(jj+1,1)=log(det(SIGMA))+np*2*log(log(t-pmax))/(t-pmax);  % HQC value
    siccrit(jj+1,1)=log(det(SIGMA))+np*log(t-pmax)/(t-pmax);         % SIC value
	
end;

% Display criterion values to verify Table 2.1
%[siccrit hqccrit aiccrit]

% Rank models for p = 0,1,2,...,pmax
[critmin,critorder]=min(aiccrit);
aichat=critorder-1;

[critmin,critorder]=min(hqccrit);
hqchat=critorder-1;

[critmin,critorder]=min(siccrit);
sichat=critorder-1;

