function DisplayEstimation(DetermTrend,Ahat,SigmahatU,p)
K = size(Ahat,1);
% Display estimation results
if DetermTrend == 0
    estimtable = table([]);
elseif DetermTrend == 1
    nuhat = Ahat(:,1);
    estimtable = table(nuhat);
elseif DetermTrend == 2
    timehat = Ahat(:,1);
    estimtable = table(timehat);
elseif DetermTrend == 3
    nuhat = Ahat(:,1);
    timehat = Ahat(:,2);
    estimtable = table(nuhat,timehat);
end
ntrend = size(estimtable,2);
Ai = reshape(Ahat(:,(1+ntrend):end),[p,K,K]);
for i = 1:p
    estimtable = [estimtable table(Ai(:,:,i),'VariableNames', {sprintf('Ahat%d',i)})];
end
estimtable = [estimtable table(SigmahatU)];
disp(estimtable)

