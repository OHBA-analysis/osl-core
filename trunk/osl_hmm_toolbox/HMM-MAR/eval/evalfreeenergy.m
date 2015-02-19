function [FrEn] = evalfreeenergy (X,T,Gamma,Xi,hmm,residuals)
% Computes the Free Energy of an HMM depending on observation model
%
% INPUT
% X            observations
% T            length of series
% Gamma        probability of states conditioned on data
% Xi           joint probability of past and future states conditioned on data
% hmm          data structure
% residuals    in case we train on residuals, the value of those.
%
% OUTPUT
% FrEn         value of the variational free energy, separated in the
%               different terms
%
% Author: Diego Vidaurre, OHBA, University of Oxford

[orders,order] = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);

Tres = sum(T) - length(T)*order;
ndim = size(X,2); K=hmm.K;
Sind = hmm.train.Sind==1; 
if ~hmm.train.zeromean, Sind = [true(1,size(Sind,2)); Sind]; end
S = hmm.train.S==1;
regressed = sum(S,1)>0;

Gammasum=sum(Gamma,1);

% compute entropy of hidden states
Entr=0;
for in=1:length(T);
    j = sum(T(1:in-1)) - order*(in-1) + 1;
    logGamma = log(Gamma(j,:));
    logGamma(isinf(-logGamma)) = log(eps);
    Entr = Entr - sum(Gamma(j,:).*logGamma);
    jj = j+1:j+T(in)-1-order;
    Gammajj = Gamma(jj,:);
    logGamma = log(Gammajj);
    logGamma(isinf(-logGamma)) = log(eps);
    Entr = Entr + sum(Gammajj(:).*logGamma(:));
end
sXi = Xi(:); 
logsXi = log(sXi);
logsXi(isinf(-logsXi)) = log(eps);
Entr = Entr - sum(sXi .* logsXi);

% Xi(Xi==0)=eps;				% avoid log(0)
% Psi=zeros(size(Xi));			% P(S_t|S_t-1)
% for k=1:K,
%     sXi=sum(squeeze(Xi(:,:,k)),2);
%     Psi(:,:,k)=Xi(:,:,k)./repmat(sXi,1,K);
% end;
% Psi(Psi==0)=eps;				% avoid log(0)
% Entr=Entr+sum(Xi(:).*log(Psi(:)),1);
%Entr=Entr+sum(Xi(:).*log(Psi(:) ./  Xi(:)),1);	% entropy of hidden states

% Free energy terms for model not including obs. model
% avLL for hidden state parameters and KL-divergence
KLdiv=dirichlet_kl(hmm.Dir_alpha,hmm.prior.Dir_alpha);
avLL = -length(T) * psi(sum(hmm.Dir_alpha));
jj = zeros(length(T),1);
for in=1:length(T);
    jj(in) = sum(T(1:in-1)) - order*(in-1) + 1;
end
for l=1:K,
    % KL-divergence for transition prob
    KLdiv=[KLdiv dirichlet_kl(hmm.Dir2d_alpha(l,:),hmm.prior.Dir2d_alpha(l,:))];
    % avLL initial state  
    avLL = avLL + sum(Gamma(jj,l)) * psi(hmm.Dir_alpha(l));
end    
% avLL remaining states  
for k=1:K,
    sXi = Xi(:,:,k); sXi = sum(sXi(:));  
    avLL = avLL - sXi * psi(sum(hmm.Dir2d_alpha(:,k)));
    for l=1:K,
        avLL = avLL + sum(Xi(:,l,k)) * psi(hmm.Dir2d_alpha(l,k));
    end;
end;

XX = formautoregr(X,T,orders,order,hmm.train.zeromean);
ltpi = ndim/2*log(2*pi); % - ndim/2;

if strcmp(hmm.train.covtype,'uniquediag')
    OmegaKL = 0;
    for n=1:ndim
        if ~regressed(n), continue; end
        OmegaKL = OmegaKL + gamma_kl(hmm.Omega.Gam_shape,hmm.prior.Omega.Gam_shape, ...
            hmm.Omega.Gam_rate(n),hmm.prior.Omega.Gam_rate(n));
    end;
    KLdiv=[KLdiv OmegaKL];
    ldetWishB=0;
    PsiWish_alphasum=0;
    for n=1:ndim,
        if ~regressed(n), continue; end
        ldetWishB=ldetWishB+0.5*log(hmm.Omega.Gam_rate(n));
        PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hmm.Omega.Gam_shape);
    end;
    C = hmm.Omega.Gam_shape ./ hmm.Omega.Gam_rate;
    avLL=avLL+ Tres*(-ltpi-ldetWishB+PsiWish_alphasum);
elseif strcmp(hmm.train.covtype,'uniquefull')
    OmegaKL = wishart_kl(hmm.Omega.Gam_rate(regressed,regressed),hmm.prior.Omega.Gam_rate(regressed,regressed), ...
        hmm.Omega.Gam_shape,hmm.prior.Omega.Gam_shape);
    KLdiv=[KLdiv OmegaKL];
    ldetWishB=0.5*logdet(hmm.Omega.Gam_rate(regressed,regressed));
    PsiWish_alphasum=0;
    for n=1:sum(regressed),
        PsiWish_alphasum=PsiWish_alphasum+psi(hmm.Omega.Gam_shape/2+0.5-n/2);   
    end;
    PsiWish_alphasum=PsiWish_alphasum*0.5;
    C = hmm.Omega.Gam_shape * hmm.Omega.Gam_irate;
    avLL=avLL+ Tres*(-ltpi-ldetWishB+PsiWish_alphasum);
end

OmegaKL = 0;

for k=1:K,
    hs=hmm.state(k);		% for ease of referencing
    pr=hmm.state(k).prior;
    
    switch hmm.train.covtype,
        case 'diag'
            ldetWishB=0;
            PsiWish_alphasum=0;
            for n=1:ndim,
                if ~regressed(n), continue; end 
                ldetWishB=ldetWishB+0.5*log(hs.Omega.Gam_rate(n));
                PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hs.Omega.Gam_shape);
            end;
            C = hs.Omega.Gam_shape ./ hs.Omega.Gam_rate;
            avLL=avLL+ Gammasum(k)*(-ltpi-ldetWishB+PsiWish_alphasum);
        case 'full'
            ldetWishB=0.5*logdet(hs.Omega.Gam_rate(regressed,regressed));
            PsiWish_alphasum=0;
            for n=1:sum(regressed),
                PsiWish_alphasum=PsiWish_alphasum+psi(hs.Omega.Gam_shape/2+0.5-n/2);   
            end;
            PsiWish_alphasum=PsiWish_alphasum*0.5;
            C = hs.Omega.Gam_shape * hs.Omega.Gam_irate;
            avLL =avLL + Gammasum(k) * (-ltpi-ldetWishB+PsiWish_alphasum);
    end;
    
    meand = zeros(size(XX,1),sum(regressed));
    if order>0
        meand = XX * hs.W.Mu_W(:,regressed);
    end
    d = residuals(:,regressed) - meand;
    
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
        Cd =  repmat(C(regressed)',1,Tres) .* d';
    else
        Cd = C(regressed,regressed) * d';
    end
    
    dist=zeros(Tres,1);
    for n=1:sum(regressed),
        dist=dist-0.5*d(:,n).*Cd(n,:)';
    end
    
    NormWishtrace=zeros(Tres,1);
    if order>0 || ~hmm.train.zeromean
        switch hmm.train.covtype,
            case {'diag','uniquediag'}
                for n=1:ndim,
                    if ~regressed(n), continue; end
                    NormWishtrace = NormWishtrace + 0.5 * C(n) *  ...
                        sum( (XX(:,Sind(:,n)) * permute(hs.W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1])) .* XX(:,Sind(:,n)), 2);
                end;
                
            case {'full','uniquefull'}
                for n1=1:ndim
                    for n2=1:ndim
                        if ~regressed(n1) || ~regressed(n2), continue; end
                        index1 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n2;
                        index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));
                        NormWishtrace = NormWishtrace + 0.5 * C(n1,n2) *  ...
                            sum( (XX(:,Sind(:,n1)) * hs.W.S_W(index1,index2)) .* XX(:,Sind(:,n2)),2);
                    end
                end
        end
    end
    
    avLL = avLL + sum(Gamma(:,k).*(dist - NormWishtrace));
    
    WKL = 0;
    if order>0 || ~hmm.train.zeromean
        if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
            for n=1:ndim
                if ~regressed(n), continue; end
                prior_prec = [];
                if hmm.train.zeromean==0
                    prior_prec = hs.prior.Mean.iS(n);
                end
                if order>0
                    ndim_n = sum(S(:,n));
                    alphamat = repmat( (hs.alpha.Gam_shape ./  hs.alpha.Gam_rate), ndim_n, 1);
                    prior_prec = [prior_prec; repmat(hs.sigma.Gam_shape(S(:,n)==1,n) ./ hs.sigma.Gam_rate(S(:,n)==1,n), length(orders), 1) .* alphamat(:)] ;
                end
                prior_prec = diag(prior_prec);
                prior_var = inv(prior_prec);
                WKL = WKL + gauss_kl(hs.W.Mu_W(Sind(:,n),n),zeros(sum(Sind(:,n)),1), ...
                    permute(hs.W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]), prior_var);
            end;
        else % full or uniquefull
            prior_prec = [];
            if hmm.train.zeromean==0
                prior_prec = hs.prior.Mean.iS;
            end
            if order>0
                sigmaterm = (hmm.state(k).sigma.Gam_shape(S==1) ./ hmm.state(k).sigma.Gam_rate(S==1) );
                sigmaterm = repmat(sigmaterm, length(orders), 1);
                alphaterm = repmat( (hs.alpha.Gam_shape ./  hs.alpha.Gam_rate), sum(S(:)), 1);
                alphaterm = alphaterm(:);
                prior_prec = [prior_prec; alphaterm .* sigmaterm];
            end
            prior_prec = diag(prior_prec);
            prior_var = inv(prior_prec);
            mu_w = hs.W.Mu_W';
            mu_w = mu_w(Sind);
            WKL = gauss_kl(mu_w,zeros(length(mu_w),1), hs.W.S_W, prior_var);
        end
    end
    
    switch hmm.train.covtype
        case 'diag'
            OmegaKL = 0;
            for n=1:ndim
                if ~regressed(n), continue; end 
                OmegaKL = OmegaKL + gamma_kl(hs.Omega.Gam_shape,hs.prior.Omega.Gam_shape, ...
                    hs.Omega.Gam_rate(n),hs.prior.Omega.Gam_rate(n));
            end;
        case 'full'
            OmegaKL = wishart_kl(hs.Omega.Gam_rate(regressed,regressed),hs.prior.Omega.Gam_rate(regressed,regressed), ...
                hs.Omega.Gam_shape,hs.prior.Omega.Gam_shape);
    end
    
    sigmaKL = 0;
    if order>0 
        for n1=1:ndim
            for n2=1:ndim
                if (hmm.train.symmetricprior && n2<n1) || S(n1,n2)==0, continue; end
                sigmaKL = sigmaKL + gamma_kl(hs.sigma.Gam_shape(n1,n2),pr.sigma.Gam_shape(n1,n2), ...
                    hs.sigma.Gam_rate(n1,n2),pr.sigma.Gam_rate(n1,n2));
            end;
        end;
    end
    
    alphaKL = 0;
    if order>0
        for i=1:length(orders)
            alphaKL = alphaKL + gamma_kl(hs.alpha.Gam_shape,pr.alpha.Gam_shape, ...
                hs.alpha.Gam_rate(i),pr.alpha.Gam_rate(i));
        end
    end
    
    KLdiv=[KLdiv OmegaKL sigmaKL alphaKL WKL];
    
    if any(isnan([OmegaKL sigmaKL alphaKL WKL])), keyboard; end
    
end;

% [errY]=hmmerror(X,T,hmm,Gamma,ones(sum(T),1),residuals);
% fprintf('XXX: %f %f %f %g \n',Entr, -avLL, +sum(KLdiv),errY(1)+errY(2))
FrEn=[Entr -avLL +KLdiv];
