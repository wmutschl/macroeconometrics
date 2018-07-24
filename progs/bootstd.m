function IRFrse = bootstd(VAR,opt)
B = 2000;
nvar = VAR.nvar;
nsteps = opt.nsteps;
const = VAR.const;
nlag = VAR.nlag;
IRFrmat=zeros(nvar,nvar,nsteps+1,B);

parfor b=1:B
    ystar = BootstrapGDP(VAR);
    VARstar = VARReducedForm(ystar,nlag,const,0);
    B0invstar = chol(VARstar.SIGu,'lower');
    IRFstar = IRFs(VARstar.Acomp,B0invstar,opt);
    % Store results
    IRFrmat(:,:,:,b)=IRFstar;
end
% Compute standard deviation across r
IRFrse=std(IRFrmat,0,4);