function ACF
% Compute acf function
R = zeros(K,K,nacf);
Gamma = zeros(K,K,nacf);   
Acomp = CompanionForm(p,Ahat);
SigmaU = zeros(K*p,K*p);
SigmaU(1:K,1:K) = SigmahatU;
GAMMA0 = reshape(inv(eye((p*K)^2)-kron(Acomp,Acomp))*SigmaU(:),[K*p K*p]);
Ai = reshape(Ahat(:,2:end),[K K p]);
Dinv=inv(diag(sqrt(diag(GAMMA0(1:K,1:K)))));
for h = 1:p
    Gamma(:,:,h) = GAMMA0(1:K,(1:K)+(h-1)*K);
    R(:,:,h)=Dinv*Gamma(:,:,h)*Dinv;
end
for h = (p+1):nacf
    for i = 1:p
        Gamma(:,:,h) = Gamma(:,:,h) + Ai(:,:,i)*Gamma(:,:,h-i);
    end
    R(:,:,h)=Dinv*Gamma(:,:,h)*Dinv;
end
% Plot autocorrelation function 

for i=1:K
    for j=1:K        
        acf = permute(R(i,j,:),[3 1 2]);
        hor=numel(acf);
        f=figure('name','Autocorrelation Theory');
        h=bar(acf,'k');
        set(gca,'xtick',1:hor);
        axis tight
        set(gca,'ylim',[-1 1]);
    end;
end;
