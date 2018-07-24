' Demonstrates how to code a 1 period rolling VAR forecast.

!n_step_forecast = 9  'Number of periods in the forecast.

'Storage, temp variables for the results

matrix(!n_step_forecast, 2) RMSE_GAP = 0
matrix(!n_step_forecast, 2) RMSE_INFL= 0
matrix(!n_step_forecast, 2) RMSE_FF = 0

scalar mse_temp_gap = 0
scalar mse_temp_infl =0
scalar mse_temp_ff = 0

for !i=0 to !n_step_forecast-1  'Main loop 

	' run a VAR on a rolling sample and create in a model object, named mod

	smpl 1980q1 1997q4+!i

     	var ch.ls 1 2 gap infl ff @ c

	delete(noerr) varmod

	ch.makemodel(varmod)

     ' use the model to run a out-of-sample forecast

	smpl 1997q4+!i 1998q1+!i 'Always 1 period ahead
	
	varmod.solveopt(s=d, d=d)
	varmod.solve

	' the difference between forecast and observation values

	delete(noerr) gap_e infl_e ff_e

	smpl 1997q4+!i 1998q1+!i

	genr gap_e = gap - gap_0
	genr infl_e = infl - infl_0
	genr ff_e = ff - ff_0

	smpl 1998q1+!i 1998q1+!i 

	group g1 gap_e infl_e ff_e

	matrix mat_diffs = g1

     	MSE_TEMP_gap = MSE_TEMP_gap + mat_diffs(1,1)^2
     	MSE_TEMP_infl = MSE_TEMP_infl + mat_diffs(1,2)^2
     	MSE_TEMP_ff = MSE_TEMP_ff + mat_diffs(1,3)^2
   
	delete mat_diffs* g1

	RMSE_gap(!i+1,1)=!i+1
	RMSE_infl(!i+1,1)=!i+1
	RMSE_ff(!i+1,1)=!i+1

	RMSE_gap(!i+1,2)=@sqrt(MSE_TEMP_gap/(!i+1))
	RMSE_infl(!i+1,2)=@sqrt(MSE_TEMP_infl/(!i+1))
	RMSE_ff(!i+1,2)=@sqrt(MSE_TEMP_ff/(!i+1))

next

show RMSE_gap
show RMSE_infl
show RMSE_ff

