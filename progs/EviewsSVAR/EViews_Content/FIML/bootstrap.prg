%s = "gali_sys" 'Estimated system; will be re-estimated to ensure that the coefficient estimates are current.
!iters = 100
!size =  0.05 'Confidence band = 1-(size/2)
!n_equations = {%s}.@neqn
delete(noerr) resid* {%s}_model
{%s}.ls
{%s}.results
{%s}.makeresids
{%s}.makemodel({%s}_model)
{%s}_model.spec
{%s}_model.addassign(i) @stochastic
smpl @all
{%s}_model.addinit(v=z) @stochastic

