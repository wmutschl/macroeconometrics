'Clean up for previous runs (if any)

'Uses chomoreno workfile.

delete(noerr) gra*

smpl 1980q1 1997q4
var var1.ls 1 2 GAP INFL FF @ c

smpl 1998q1 2000q1

'Use the VAR built-in forecasting routine to produce a nine period ahead dynamic forecast (default).

var1.forecast _f  'Forecasts will be contained in gap_f, infl_f and ff_f.

show gap gap_f infl infl_f ff ff_f

Table summary_results

summary_results.setwidth(1) 35

summary_results(1,3) = "RMSE"
summary_results(1,4) = "MAE"
summary_results(1,5) = "MAPE"
summary_results(1,6) = "THEIL"
summary_results(1,7) = "THEIL_2"

summary_results(2,1) = "Model 1: Unrestricted"
summary_results(2,2) = "infl"
summary_results(2,3) =  @rmse(infl,infl_f)
summary_results(2,4) =  @mae(infl,infl_f)
summary_results(2,5) =  @mape(infl,infl_f)
summary_results(2,6) = @theil(infl,infl_f)
summary_results(2,7) = @theil2(infl,infl_f)

summary_results(3,1) = "Model 1: Unrestricted"
summary_results(3,2) = "gap"
summary_results(3,3) = @rmse(gap,gap_f)
summary_results(3,4) = @mae(gap,gap_f)
summary_results(3,5) =  @mape(gap,gap_f)
summary_results(3,6) = @theil(gap,gap_f)
summary_results(3,7) = @theil2(gap,gap_f)

summary_results(4,1) = "Model 1: Unrestricted"
summary_results(4,2) = "ff"
summary_results(4,3) = @rmse(ff,ff_f)
summary_results(4,4) = @mae(ff,ff_f)
summary_results(4,5) =  @mape(ff,ff_f)
summary_results(4,6) = @theil(ff,ff_f)
summary_results(4,7) = @theil2(ff,ff_f)

show summary_results 

'Create a model from the estimated VAR

delete(noerr) Model_1

var1.makemodel(Model_1)

show Model_1

'change sample to forecast period
smpl 1998q1 2000q1

'solve model to obtain dynamic forecasts
Model_1.solve  'Solution will be in _0 variables....

show gap gap_0 infl infl_0 ff ff_0

show (gap_f - gap_0) (infl_f-infl_0) (ff_f - ff_0) 'This should be zero; the 2 ways of generating dynamic forecasts are the same.

Table summary_results_o

summary_results.setwidth(1) 35

summary_results_o(1,3) = "RMSE"
summary_results_o(1,4) = "MAE"
summary_results_o(1,5) = "MAPE"
summary_results_o(1,6) = "THEIL"
summary_results_o(1,7) = "THEIL_2"

summary_results_o(2,1) = "Model 1: Unrestricted"
summary_results_o(2,2) = "infl"
summary_results_o(2,3) = @rmse(infl,infl_0)
summary_results_o(2,4) = @mae(infl,infl_0)
summary_results_o(2,5) = @mape(infl,infl_0)
summary_results_o(2,6) = @theil(infl,infl_0)
summary_results_o(2,7) = @theil2(infl,infl_0)

summary_results_o(3,1) = "Model 1: Unrestricted"
summary_results_o(3,2) = "gap"
summary_results_o(3,3) = @rmse(gap,gap_0)
summary_results_o(3,4) = @mae(gap,gap_0)
summary_results_o(3,5) = @mape(gap,gap_0)
summary_results_o(3,6) = @theil(gap,gap_0)
summary_results_o(3,7) = @theil2(gap,gap_0)

summary_results_o(4,1) = "Model 1: Unrestricted"
summary_results_o(4,2) = "ff"
summary_results_o(4,3) = @rmse(ff,ff_0)
summary_results_o(4,4) = @mae(ff,ff_0)
summary_results_o(4,5) = @mape(ff,ff_0)
summary_results_o(4,6) = @theil(ff,ff_0)
summary_results_o(4,7) = @theil2(ff,ff_0)

show summary_results_o

'Now let's do a conditional forecast...by conditioning on the federal funds rate.  There is simpler way of doing this by dropping the FF equation from the model and adding FF=FF_1 (see forecast_conditional for an example).

Model_1.scenario(c) "actuals"
Model_1.scenario "scenario 1" 'Must use an alternative scenario to calculate a conditional forecast
'mod1.override ff 'This will use ff_1
'mod1.solve
'mod1.vars
Model_1.exclude(m) ff
Model_1.override(m) ff

smpl @all
 
genr ff_1=ff_cond

show ff_1

smpl 1998q1 2000q1

Model_1.solve

show gap_0 gap_1 (gap_0-gap_1)  'Third variable shows the difference.

' plot actual and forecasts

smpl 1996:1 2000:1

   group gtmp  GAP gap_1
   freeze(gra1) gtmp.line
   %gname = %gname + "gra" + @str(1) + " "

 ' plot actual and forecasts
smpl 1996:1 2000:1

   group gtmp  INFL INFL_1
   freeze(gra2) gtmp.line
   %gname = %gname + "gra" + @str(2) + " "

 ' plot actual and forecasts
smpl 1996:1 2000:1

   group gtmp  FF FF_1
   freeze(gra3) gtmp.line
   %gname = %gname + "gra" + @str(3) + " "

' merge all graphs into one
freeze(graph_forcast) {%gname}
graph_forcast.options size(8,2)
graph_forcast.align(1, 0.1, 0.5)
graph_forcast.legend position(0.1,0.1)
'gfcst.scale(left) +zeroline
graph_forcast.draw( dashline,left,rgb(155,155,155) ) 0.0
 
show graph_forcast
 smpl 1998q1 2000q1
show graph_forcast

Table summary_results_c

summary_results.setwidth(1) 35

summary_results_c(1,3) = "RMSE"
summary_results_c(1,4) = "MAE"
summary_results_c(1,5) = "MAPE"
summary_results_c(1,6) = "THEIL"
summary_results_c(1,7) = "THEIL_2"

summary_results_c(2,1) = "Model 1: Conditional"
summary_results_c(2,2) = "infl"
summary_results_c(2,3) = @rmse(infl,infl_1)
summary_results_c(2,4) = @mae(infl,infl_1)
summary_results_c(2,5) = @mape(infl,infl_1)
summary_results_c(2,6) = @theil(infl,infl_1)
summary_results_c(2,7) = @theil2(infl,infl_1)

summary_results_c(3,1) = "Model 1: Conditional"
summary_results_c(3,2) = "gap"
summary_results_c(3,3) = @rmse(gap,gap_1)
summary_results_c(3,4) = @mae(gap,gap_1)
summary_results_c(3,5) = @mape(gap,gap_1)
summary_results_c(3,6) = @theil(gap,gap_1)
summary_results_c(3,7) = @theil2(gap,gap_1)

summary_results_c(4,1) = "Model 1: Conditional"
summary_results_c(4,2) = "ff"
summary_results_c(4,3) = @rmse(ff,ff_1)
summary_results_c(4,4) = @mae(ff,ff_1)
summary_results_c(4,5) = @mape(ff,ff_1)
summary_results_c(4,6) = @theil(ff,ff_1)
summary_results_c(4,7) = @theil2(ff,ff_1)

show summary_results_c


