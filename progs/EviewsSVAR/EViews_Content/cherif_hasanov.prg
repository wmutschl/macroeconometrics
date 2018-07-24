'uses debt.wf1
'Replicates the results presented in Section 4.6.4
smpl 1981q3 2007q4
b_var_alt.ls 1 2 pb dy dp in dhat  @ c debt(-1) debt(-2) 
b_var_alt.svar(rtype=patsr, aname=a_mat_alt, bname=b_mat_alt, f0=u)
b_var_alt.impulse(100, imp=struct, se=a) dhat @ 1 2 3 4 5


