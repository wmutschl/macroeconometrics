function [StructMod,opt] = SVECM(VECM,opt)
% =======================================================================
% Perform and display summary of Augmented Dickey Fuller test for variables
% in levels as well as in first differences
% =======================================================================
% SVECM(VECM,opt)
% -----------------------------------------------------------------------
% INPUT
% -----------------------------------------------------------------------
% OUTPUT
% -----------------------------------------------------------------------
% CALL
% =======================================================================
% Willi Mutschler, March 2017
% willi@mutschler.eu
algorithm = opt.algorithm; % String for algorithm used for estimation, values are 'RWZ', 'Scoring', 'Optimization'
StartValueMethod = opt.StartValueMethod;

nlag = opt.nlag;
K = opt.nvars;
T = size(VECM.res,1);
SIGMAUHAT = VECM.EstCov;
alph = VECM.paramVals.A;
betta = VECM.paramVals.B;
Rshort = opt.Restr.Rshort;
Rlong = opt.Restr.RlongVECM;
Rsign = opt.Restr.Rsign;

I_K = eye(K);
I_KK = eye(K^2);
GAM = I_K;
if nlag >0
    for i = 1:(nlag-1)
        GAM = GAM - VECM.paramVals.(sprintf('B%d',i));
    end
end
betta_o = null(betta');	%/* orthogonal complement of beta */
alph_o = null(alph');	%/* orthogonal complemt of alpha  */
XI = betta_o*((alph_o'*GAM*betta_o)\transpose(alph_o)); %/* total impact matrix */

selB0inv   = find(isnan(Rshort)==0);
selUpsilon = find(isnan(Rlong)==0);

switch algorithm
    case 'Cholesky'
        B0inv = chol(SIGMAUHAT,'lower');
    case 'Optimization'
        % Numerical Optimization
        TolX = 1e-4;
        TolFun = 1e-9;
        MaxFunEvals = 1e+50000;
        MaxIter = 1000;
        NumIter = 3;
        OptimAlgorithm = 'trust-region-dogleg';
        %Display = 'off';        
        if StartValueMethod == 1
            B0inv_0= SIGMAUHAT^.5; 
        elseif StartValueMethod == 2
            B0inv_0 = chol(SIGMAUHAT,'lower');
        end
        options=optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'MaxIter',MaxIter,'Algorithm',OptimAlgorithm);
        [B0inv,fval,exitflag,output]=fsolve('f_SVECM',B0inv_0,options,SIGMAUHAT,XI,Rshort,Rlong,selB0inv,selUpsilon);
        if strcmp(output.algorithm,'levenberg-marquardt')
            options.Algorithm = output.algorithm;            
        end
        for r = 1:NumIter
            [B0inv,fval,exitflag,output]=fsolve('f_SVAR',B0inv,options,SIGMAUHAT,XI,Rshort,Rlong,selB0inv,selUpsilon);
        end

    case 'Scoring'
        % Set additional options
        eps1_tol = 1E-6;	   %tolerance of relative change in parameters 
        eps2_tol = 1E-10; 	   %tolerance for relative change in logLik
        maxIterations = 500;   %maximum number of iterations
        StartValueMethod = 1; %1: draw randomly, 2: fixed values set by "fixStart"
        imaxRetries	= 10;      %max. retries of randomly drawing starting values 
        fixStart = .1;         %fixed starting value, only needed when iStartValueMethod  = 2         
        iCorr = 0;             %decompose correlation matrix first to obtain SV
       
        R_B0inv = I_KK(selB0inv,:);
        R_Upsilon = I_KK(selUpsilon,:);

        if isempty(R_Upsilon);
            R_Upsilon = [];
        end	
        R_Upsilon = R_Upsilon*kron(I_K,XI);

        if isempty(R_B0inv) == 0
            R = [R_Upsilon;R_B0inv];   	
        else
            R = R_Upsilon;
        end     	
   	
        Sb = null(R);
        S = [zeros(K^2,size(Sb,2));Sb];
        vs = [vec(eye(K));zeros(K^2,1)];
        nfreeall = size(Sb,2);

        if isempty(R_Upsilon) == 0
            nResUpsilon  = rank(R_Upsilon);
        else
            nResUpsilon = 0;
        end

        if isempty(R_B0inv) == 0;
            nResB0inv = rank(R_B0inv);
        else
            nResB0inv = 0;
        end

        %Starting Values for algorithm
        if StartValueMethod == 1
            vg = randn(size(S,2),1)*.1;
        elseif StartValueMethod == 2
            vg = ones(size(S,2),1)*fixStart;
        end
        
        ident = 0;
        if (nResB0inv+nResUpsilon) < (K*(K-1)/2)
            DegreeofOverIdent = 0;
        else 	 
            vAB = S*vg+vs;
            A = reshape(vAB(1:K^2),K,K)';
            B = reshape(vAB(K^2+1:2*K^2),K,K)';
            Com = commutation(K,K);
            Binv = B\I_K;
            V   = [(I_KK+Com)*kron(((A\B))',Binv)  -1*(I_KK+Com)*kron(I_K,Binv)];
            DegreeofOverIdent = K*(K+1)/2 - nfreeall;

            if isempty(V)==0 && isempty(S)==0
                if sum(eig((V*S)'*(V*S)) < 1E-10) == 0
                    ident = 1;
                end
            end	  	
        end

        if ident == 1;     %/* only proceed when model is identified, otherwise inform user*/    		
            if iCorr == 1;	 %/* decompose correlation matrix first to obtain SV 	*/
                [vg_corr,i_iRT,i_itercorr, i_coCorr]  = DecompCorrelation_svar_var(mSigmaU,mS,vs,vg,T,eps1_tol,eps2_tol,DegreeofOverIdent,maxIterations,model, imaxRetries);	
                vg = vg_corr; 
                if i_iRT == 0 && i_coCorr == 0;	%/* only estimate when 1st step was successful*/  
                    [mB,vg_b,mSigmaUt,i,LLnew,LR_stat,LR_prob, noConvergence] = MLEstSVECM_svec_var(SIGMAUHAT, S,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,3);
                end
            else
                %[B0inv,vg_b,SigmaUt,i,LLnew,LR_stat,LR_prob, noConvergence] =   MLEstSVECM_svec_var(SIGMAUHAT, S,vs,vg,T,eps1_tol, eps2_tol, DegreeofOverIdent,maxIterations,3);
                maxls   = 1;
                vecAB   = S*vg+vs;
                A   	= reshape(vecAB(1:K^2,1),K,K);
                B   	= reshape(vecAB(K^2+1:2*K^2,1),K,K);
                Com 	= commutation(K,K);
                noConvergence = 0;

                vgold   = vg;
                eps1 	= 100;
                eps2 	= 100;

                %mSigmaUt = inv(A)*B*B'*inv(A');
                AinvB=A\B;
                SigmaUt = AinvB*AinvB';
                mK = B\A;
                llold = T/2*log(det(mK)^2)-T/2*sum(diag((mK'*mK*SigmaUt)));

                i = 0;
                while (isempty(find(eps1 > eps1_tol)) || (eps2 > eps2_tol)) && i<maxIterations		
                    mK = B\A;            
                    IAB = T*[ kron(inv(mK),inv(B')) ; -1*kron(I_K,inv(B')) ]*(I_KK+ Com)*[ kron(inv(mK'),inv(B)) -1*kron(I_K,inv(B))];
                    Iga = S'*IAB*S;
                    v_scoreK 	 = T*(vec(inv(mK)')'-vec(mK)'*kron(SIGMAUHAT,I_K));
                    v_scoreAB 	 = v_scoreK * [kron(I_K,inv(B))  -1*kron(A'*inv(B'),inv(B))];
                    v_scoregamma = v_scoreAB*S;  

                   tmp = inv(Iga);
                   length = max(abs(tmp*v_scoregamma'));
                   if length > maxls
                    lambda = maxls/length;
                   else
                     lambda =1;
                   end   

                   vg = vgold +lambda*tmp*v_scoregamma';

                    vecAB = S*vg+vs;

                    A = reshape(vecAB(1:K^2),K,K);
                    B = reshape(vecAB(K^2+1:2*K^2),K,K);   
                    mK = inv(B)*A;
                    SigmaUt = inv(A)*B*B'*inv(A');

                    llnew  = T/2*log(det(mK)^2)-T/2*sum(diag((mK'*mK*SigmaUt)));

                    eps2 = abs((llnew-llold)/llold);
                    eps1 = abs((vg-vgold)./vgold);
                    vgold = vg;	
                    llold  = llnew;
                    i = i +1;  
                end
                B0inv =inv(A)*B;

                if isempty(find(eps1>eps1_tol))==0 || (eps2>eps2_tol);
                  noConvergence = 1;
                end

                vg_B = inv(S'*S)*S'*vec([A B]);

                %/* compute LR test for over-identifying restr. */				    
                if DegreeofOverIdent > 0;       
                    LR_stat =  T*(log(det(SigmaUt))-log(det(SIGMAUHAT)));
                    LR_prob =  chi2cdf(LR_stat,DegreeofOverIdent);
                else
                    LR_stat = 0;
                    LR_prob = 0;
                end
            end
        end
        
        %vg_BS 	  = vg_b;

    %     %/* Std. deviations and t-ratios */
    % 
    %     m_se_B=[];
    %     m_tv_B=[];
    %     m_se_mCB=[];
    %     m_tv_mCB=[];
    % 
    %     %/* Bootstrap std deviations */
    %     if BS_int  
    %       [m_se_B,m_se_mCB,m_tv_B,m_tv_mCB] = BootstrapStdErr_svec_var(var,beta,beta_d,bootRep,seed,mB_Res,mC1_Res,vg_bs,T,K,mB_point,mCB_point,eps1_tol, eps2_tol,DegreeofOverIdent, maxIterations);	
    %     end

    case 'RWZ'
        if StartValueMethod == 1
            B0inv_0= SIGMAUHAT^.5; 
        elseif StartValueMethod == 2
            B0inv_0 = chol(SIGMAUHAT,'lower');
        end
        
        f = [Rshort;
             Rlong];
        f(isnan(f)) = 1;

        [Q,index,flag] = findQs(K,f);
%        TR = 500; draws = 1000;
%        shocks = eye(K);
%        shock_pos = logical(shocks); % position of shock
%        R = zeros(K,TR,draws,K); % Contains impulse resonse functions
        % 1st dimension = variable
        % 2nd dimension = time
        % 3rd dimension = draw
        % 4th dimension = shock
        draws = 1;
        counter = 1;

        while counter < draws+1
            B0inv = generateDraw(B0inv_0,K);
            %P = findP(C,B,Q,p,k,index);
            L0 = B0inv;
            LUpsilon = XI*B0inv;
            F = [L0;LUpsilon];
            P = zeros(K,K);    
            for ii = 1:K
                if ii == 1
                    Qtilde = Q{ii}*F;
                else
                    Qtilde = [Q{ii}*F;P'];
                end
                [QQ,RR] = qr(Qtilde');
                P_temp = QQ(:,end);
                P(:,ii) = P_temp;    
            end
            P = P(:,index);
%            W = B0inv*P;
% 
%             for jj = 1:length(shocks)
% 
%                 if W(var_pos(jj),jj) < 0
%                     shock = -shocks(:,jj);
%                 else
%                     shock = shocks(:,jj);
%                 end
% 
%                 V = zeros(k*p,T);
% 
%                 V(1:k,1) = W*shock;
% 
%                 chk = W*shock;
%                 sr_index = ~isnan(sr(:,jj));
%                 tmp = sign(chk(sr_index)) - sr(sr_index,jj);
% 
%                 if any(tmp~=0)
%                     jj = 0;
%                     break
%                 end
% 
%                 for ii = 2:T
%                     V(:,ii) = alpha*V(:,ii-1);
%                 end
% 
%                 R(:,:,counter,jj) = V(1:k,:);
% 
%             end
% 
%             if jj == length(shocks)
                counter = counter + 1;
%            end
        end
        B0inv = L0*P(:,index);
end

if sum(diag(B0inv)<0) ~= 0               %/* normalize sign of mB */
	x = diag(B0inv)<0;
    B0inv(:,find(x==1)) = -1*B0inv(:,find(x==1));
end

StructMod = VECM;
StructMod.B0inv_point = B0inv;
StructMod.Upsilon_point = XI*B0inv;

