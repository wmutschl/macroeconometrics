‘'This command opens an Eviews Workfile with Cho and Moreno data

'wfopen "U:\My Documents\CURSOS\MFORE\OPR\BVAR\chom.WF1"

'wfopen "U:\My Documents\CURSOS\MFORE\OPR\TEXT\ChoMor\chom_forecast.WF1"

smpl 1980q1 1997q4

'This command creates a group named y with the time series we will use in the VAR:d(gap) (change in the output gap) infl (inflation) and d(ff) (change in the Fed funds interest rate)

group y d(gap) infl d(ff)

''lag exclusion: this command runs a VAR with the variables included in the group Y and a constant. It includes six lags

var chounrestr.ls 1 6 y @ c

show chounrestr

'We run a lag exclusion test to know which is the best number of lags to be included in the VAR. In addition we use the command freeze to generate a table with the results 

delete(noerr) lags tab1

freeze(tab1) chounrestr.testlags(name=lags) 

show tab1

var chounrestr.ls 1 2 y @ c

chounrestr.makemodel(chounres)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'This command estimates a Bayesian VAR using the Litterman (Minnesota) priors

'Bayesian VAR with 2 lags:

var chobvar.bvar(prior=lit,initcov=uni,df,mu1=0, L1=0.1,L2=0.99,L3=1) 1 2 d(gap) infl d(ff) @ c

chobvar.makemodel(chobv)

'MU1 is the mean of the coefficient of the first lag of the dependent variables.

'MU1 was set equal to 0 because the theory would tell us that these variables should be stationary. 

' L1 measures the overall tightness of the priors. It was set equal to 0.5 meaning that we have a strong belief in the priors. 

' If we had no clue about the priors, we could have set that number to 10, which is a looser beilief about the priors. So the priors would be less binding and.the estimation would be closer to the unrestricted VAR

' L2 indicates the importance of lags of the other variables in each equation. When L2 is very small the other lagged variables play a less important role in the ecuations, so their coefficients go to zero.

' L2 is a number between 0 and.1. 

' L3 is the decay of the importance of the lags larger than one in each equation. A higher L3 will imply a very fast decay of the coefficients of the lags.

'show chobvar.impulse(20, m) d(gap) infl d(ff) @ d(gap) infl d(ff)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'The following command estimates a BVAR using the NORMAL-WISHART PRIOR

var chobvar1.bvar(prior=nw,df,mu1=0, L1=.1) 1 2 d(gap) infl d(ff) @ c

chobvar1.makemodel(chobv1)

'EViews does not allow the setting of the prior for the variance-covarianze of the residuals. The options here are to set MU, which governs again the prior for the series' process: stationary?

' I(1)? If you choose a small MU (less than 1) your prior is that those series are stationary. If MU is set equial to 1, your prior is that the series have a unit root.

'The other hyper parameter you are allowed to choose in the case of a Normal-Wishart is the overall tightness. Here EViews is doing the opposite to what is said in the manual.

'In this case a large L1 implies that you have great certainty about your prior. Here EViews works the opposite to how the Minnesota prior works.

'show chobvar1.impulse(20, m) gap infl ff @ d(gap) infl d(ff) 

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'The following command estimates a BVAR using Sims-Zha Normal-Wishart (dummy approach).'This approach allows for the existence of cointegration among I(1) variables in the VAR.

var chobvar2.bvar(prior=sznw,df,mu5=0.00,mu6=0.00, L1=0.1, L3=1, L0=1) 1 2 d(gap) infl d(ff) @ c

chobvar2.makemodel(chobv2)

'if we increase MU5 we will shrink the coefficient of the lag towards 1, meaning that the series will have a unit root. 

'If we increase MU6 we allow for the existence of a cointegrating vector.Again, a high MU6 pushes the system towards having a cointegrating vector.

'by setting them equal to 0, we are trying to shut them down. Thus, they do not play a role and the prior will be similar to the Normal-Wishart prior.

'L1 controls the overall tightness as before.

'L3 governs the dacay of coefficients for higher lags. 

'L0 is related to the tightness of the prior on the residual Variance-Covariance matrix

'Conceptually,Sims and Zha developed this type of prior for a case in which the variables have unit roots.

'show chobvar2.impulse(20, m) d(gap) infl d(ff) @ d(gap) infl d(ff) 

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

smpl 1998.1  2000.1

chounres.solve

genr aux1 = @mean( (infl-infl_0)^2)^(0.5)

genr aux2 = @mean(@abs(infl-infl_0))

genr aux3 = @mape(infl,infl_0)

genr aux4 =(@mean( (infl-infl_0)^2)^(0.5))/(@mean( (infl)^2)^(0.5)+@mean( (infl_0)^2)^(0.5))

 
scalar  RMSE =  @last(aux1)

scalar  MAE    =   @last(aux2)

scalar MAPE  =   @last(aux3)

scalar  Theil   =   @last(aux4)

 
genr aux1a = @mean( (gap-gap_0)^2)^(0.5)

genr aux2a = @mean(@abs(gap-gap_0))

genr aux3a  = @mape(gap, gap_0)

genr aux4a =(@mean( (gap-gap_0)^2)^(0.5))/(@mean( (gap)^2)^(0.5)+@mean( (gap_0)^2)^(0.5))

 
scalar  RMSEa =  @last(aux1a)

scalar  MAEa    =   @last(aux2a)

scalar MAPEa  =   @last(aux3a)

scalar  Theila   =   @last(aux4a)

chobv.solve

genr aux1 = @mean( (infl-infl_0)^2)^(0.5)

genr aux2 = @mean(@abs(infl-infl_0))

genr aux3 = @mape(infl,infl_0)

genr aux4 =(@mean( (infl-infl_0)^2)^(0.5))/(@mean( (infl)^2)^(0.5)+@mean( (infl_0)^2)^(0.5))

scalar  RMSE1 =  @last(aux1)

scalar  MAE1    =   @last(aux2)

scalar MAPE1  =   @last(aux3)

scalar  Theil1   =   @last(aux4)

genr aux1a = @mean( (gap-gap_0)^2)^(0.5)

genr aux2a = @mean(@abs(gap-gap_0))

genr aux3a  = @mape(gap, gap_0)

genr aux4a =(@mean( (gap-gap_0)^2)^(0.5))/(@mean( (gap)^2)^(0.5)+@mean( (gap_0)^2)^(0.5))

 

scalar  RMSE1a =  @last(aux1a)

scalar  MAE1a    =   @last(aux2a)

scalar MAPE1a  =   @last(aux3a)

scalar  Theil1a   =   @last(aux4a)

 

 

chobv1.solve

genr aux1 = @mean( (infl-infl_0)^2)^(0.5)

genr aux2 = @mean(@abs(infl-infl_0))

genr aux3 = @mape(infl,infl_0)

genr aux4 =(@mean( (infl-infl_0)^2)^(0.5))/(@mean( (infl)^2)^(0.5)+@mean( (infl_0)^2)^(0.5))

 

scalar  RMSE2 =  @last(aux1)

scalar  MAE2    =   @last(aux2)

scalar MAPE2  =   @last(aux3)

scalar  Theil2   =   @last(aux4)

 

genr aux1a = @mean( (gap-gap_0)^2)^(0.5)

genr aux2a = @mean(@abs(gap-gap_0))

genr aux3a  = @mape(gap, gap_0)

genr aux4a = (@mean( (gap-gap_0)^2)^(0.5))/(@mean( (gap)^2)^(0.5)+@mean( (gap_0)^2)^(0.5))

 

scalar  RMSE2a =  @last(aux1a)

scalar  MAE2a    =   @last(aux2a)

scalar MAPE2a  =   @last(aux3a)

scalar  Theil2a  =   @last(aux4a)

 

chobv2.solve

genr aux1 = @mean( (infl-infl_0)^2)^(0.5)

genr aux2 = @mean(@abs(infl-infl_0))

genr aux3 = @mape(infl,infl_0)

genr aux4 = (@mean( (infl-infl_0)^2)^(0.5))/(@mean( (infl)^2)^(0.5)+@mean( (infl_0)^2)^(0.5))

 

scalar  RMSE3 =  @last(aux1)

scalar  MAE3    =   @last(aux2)

scalar MAPE3  =   @last(aux3)

scalar  Theil3   =   @last(aux4)

 

genr aux1a = @mean( (gap-gap_0)^2)^(0.5)

genr aux2a = @mean(@abs(gap-gap_0))

genr aux3a  = @mape(gap, gap_0)

genr aux4a = (@mean( (gap-gap_0)^2)^(0.5))/(@mean( (gap)^2)^(0.5)+@mean( (gap_0)^2)^(0.5))

 

scalar  RMSE3a  =  @last(aux1a)

scalar  MAE3a     =   @last(aux2a)

scalar  MAPE3a   =   @last(aux3a)

scalar  Theil3a    =   @last(aux4a)

 

Table summary_results

 

summary_results.setwidth(1) 35

summary_results(1,3) = "RMSE"

summary_results(1,4) = "MAE"

summary_results(1,5) = "MAPE"

summary_results(1,6) = "THEIL"

 

summary_results(2,1) = "Model 1: Unrestricted"

summary_results(2,2) = "infl"

summary_results(2,3) =  RMSE

summary_results(2,4) =  MAE

summary_results(2,5) =  MAPE

summary_results(2,6) = THEIL

 

summary_results(3,1) = "Model 1: Unrestricted"

summary_results(3,2) = "gap"

summary_results(3,3) =  RMSEa

summary_results(3,4) =  MAEa

summary_results(3,5) =  MAPEa

summary_results(3,6) = THEILa

 

summary_results(5,1) = "Model 2: Minnesota"

summary_results(5,2) = "infl"

summary_results(5,3) = RMSE1

summary_results(5,4) = MAE1

summary_results(5,5) = MAPE1

summary_results(5,6) = THEIL1

 

summary_results(6,1) = "Model 2: Minnesota"

summary_results(6,2) = "gap"

summary_results(6,3) = RMSE1a

summary_results(6,4) = MAE1a

summary_results(6,5) = MAPE1a

summary_results(6,6) = THEIL1a

 

summary_results(8,1) = "Model 3: Normal-Wishart"

summary_results(8,2) = "infl"

summary_results(8,3) = RMSE2

summary_results(8,4) = MAE2

summary_results(8,5) = MAPE2

summary_results(8,6) = THEIL2

 

summary_results(9,1) = "Model 3: Normal-Wishart"

summary_results(9,2) = "gap"

summary_results(9,3) = RMSE2a

summary_results(9,4) = MAE2a

summary_results(9,5) = MAPE2a

summary_results(9,6) = THEIL2a

 

summary_results(11,1) = "Model 4: Sims-Zha_Normal-Wishart"

summary_results(11,2) = "infl"

summary_results(11,3) = RMSE3

summary_results(11,4) = MAE3

summary_results(11,5) = MAPE3

summary_results(11,6) = THEIL3

 

summary_results(12,1) = "Model 4: Sims-Zha_Normal-Wishart"

summary_results(12,2) = "gap"

summary_results(12,3) = RMSE3a

summary_results(12,4) = MAE3a

summary_results(12,5) = MAPE3a

summary_results(12,6) = THEIL3a

 

 

show summary_results


