function mB_BS = GetStructDecompinBS_svec_var(mSigmaU, mS, vs, vg, T, eps1_tol, eps2_tol, DegreeofOverIdent, maxIterations,model,j)
% /*
% **  MLEstSVECM to be used in boostrap!
% **	
% **	- if no convergence, different starting values are used
% **	
% **	- stops if no convergence after 10 retries with different starting values
% */

noCon_BS = 1;     		          
c_count  = 1;
while noCon_BS == 1
    [mB_BS,tmp,tmp,i,tmp,tmp,tmp,noCon_BS] = MLEstSVECM_svec_var(mSigmaU,mS,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,model);		
    if noCon_BS  == 1
        fprintf('no convergence in repl. %d trying different start values...;\n',j);
        vg = rndn(size(vg,1),1);
    end;
		
    c_count = c_count + 1;
    if c_count > 10;
        mB_BS = missex(mB_BS,ones(rows(mB_BS),cols(mB_BS)));
        error('No convergence in structural decomposition after 50 retries!');        
    end
end
