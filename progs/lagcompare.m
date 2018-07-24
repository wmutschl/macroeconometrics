M = 5;

if M ==1
    % Set up model 1
    A1 = [0.5, 0; 
          0, 0.5];
    A2 = [0.1, 0; 
          0, 0.1];
    SIGu = eye(2);
    p=2;
elseif M==2
    % Set up model 2
    A1 = [0.2, 0.1; 
          -0.3, 0.5];
    A2 = [0.02, 0; 
          0, 0.05];
    SIGu = eye(2);
    p=2;
elseif M==3
    % Set up model 3
    A1 = [0.5, -0.1; 
          -0.01, 0.5];
    A2 = [-0.3, 0.1; 
          0, -0.3];
    SIGu = eye(2);
    p=2;
elseif M==4
    % Set up model 3
    A1 = [0.1, -0.1; 
          0, 0.1];
    A2 = [0.7, 0.3; 
          0.2, 0.4];
    SIGu = eye(2);
    p=2;
elseif M==5
    % Set up model 3
    A1 = [2.4, 1.0; 
          0, 1.1];
    A2 = [-2.15, -0.9; 
          0, -0.41];
    A3 = [0.852, 0.2; 
          0, 0.06];
    A4 = [-0.126, 0; 
          0, 0.0003];
    SIGu = eye(2);
    p = 4;
end
SIGu = [0.9, 0.2;
        0.2, 0.5];
T = [80 160 240 500 10000];
nT = length(T);
burnin = 100;
y = zeros(2,T(nT)+burnin);
pmax = 8;
R = 100;
sic = nan(R,nT); hqc=nan(R,nT);aic = nan(R,nT);
parfor_progress(R); % Initialize 
for r=1:R
    for t = 1+p:T(end)+burnin
        if M ==5
            y(:,t) = A1*y(:,t-1) + A2*y(:,t-2) + A3*y(:,t-3) + A4*y(:,t-4) + transpose(mvnrnd([0;0],SIGu));
        else
            y(:,t) = A1*y(:,t-1) + A2*y(:,t-2) + transpose(mvnrnd([0;0],SIGu));
        end
    end
    for j = 1:nT
        ENDO = transpose(y(:,end-T(j)+1:end));
        [sic(r,j),hqc(r,j),aic(r,j)]=pfind(ENDO,pmax);
    end
    parfor_progress; % Count 
end
parfor_progress(0); % Clean up

tabulate(aic(:,1))
tabulate(hqc(:,1))
tabulate(sic(:,1))

tabulate(aic(:,2))
tabulate(hqc(:,2))
tabulate(sic(:,2))

tabulate(aic(:,3))
tabulate(hqc(:,3))
tabulate(sic(:,3))

tabulate(aic(:,4))
tabulate(hqc(:,4))
tabulate(sic(:,4))

tabulate(aic(:,5))
tabulate(hqc(:,5))
tabulate(sic(:,5))