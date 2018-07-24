function y_artificial = SimulateData(VAR,u,opt)
%% STEP 2.1: initial values for the artificial data
% Intialize the first nlag observations with real data
nlag = opt.nlag;
nvars = opt.nvars;
nobs = VAR.nobs;
y_artificial = zeros(nobs+nlag,nvars);
LAG=[];
for jj = 1:nlag
    y_artificial(jj,:) = ENDO(jj,:);
    LAG = [y_artificial(jj,:) LAG]; 
end
% Initialize the artificial series and the LAGplus vector
T = [1:nobs]';
if const==0
    LAGplus = LAG;
elseif const==1
    LAGplus = [1 LAG];
elseif const==2
    LAGplus = [1 T(1) LAG]; 
elseif const==3
    T = [1:nobs]';
    LAGplus = [1 T(1) T(1).^2 LAG];
end
if nvars_ex~=0
    LAGplus = [LAGplus StructMod.X_EX(jj-nlag+1,:)];
end

%% STEP 2.2: generate artificial series
% From observation nlag+1 to nobs, compute the artificial data
for jj = nlag+1:nobs+nlag
    for mm = 1:nvars
        % Compute the value for time=jj
        y_artificial(jj,mm) = LAGplus * Ft(1:end,mm) + u(jj-nlag,mm);
    end
    % now update the LAG matrix
    if jj<nobs+nlag
        LAG = [y_artificial(jj,:) LAG(1,1:(nlag-1)*nvars)];
        if const==0
            LAGplus = LAG;
        elseif const==1
            LAGplus = [1 LAG];
        elseif const==2
            LAGplus = [1 T(jj-nlag+1) LAG];
        elseif const==3
            LAGplus = [1 T(jj-nlag+1) T(jj-nlag+1).^2 LAG];
        end
        if nvars_ex~=0
            LAGplus = [LAGplus StructMod.X_EX(jj-nlag+1,:)];
        end
    end
end
