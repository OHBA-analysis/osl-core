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


if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end); 
else orders = []; order = 0; end 

Tres = sum(T) - length(T)*order;
N = length(T);
ndim = size(X,2);
K=hmm.K;

Gammasum=sum(Gamma,1);

% compute entropy of hidden states
% Entropy of initial state
Entr=0;
for in=1:length(T);
    j = sum(T(1:in-1)) - order*(in-1) + 1;
    Entr = Entr + sum((Gamma(j,:)+eps).*log(Gamma(j,:)+eps),2);
end
Xi=Xi+eps;				% avoid log(0)
Psi=zeros(size(Xi));			% P(S_t|S_t-1)
for k=1:K,
    sXi=sum(squeeze(Xi(:,:,k)),2);
    Psi(:,:,k)=Xi(:,:,k)./repmat(sXi,1,K);
end;
Psi=Psi+eps;				% avoid log(0)
Entr=Entr+sum(Xi(:).*log(Psi(:)),1);	% entropy of hidden states

% Free energy terms for model not including obs. model
% avLL for hidden state parameters and KL-divergence
avLL=0;
PsiDir2d_alpha=zeros(K,K);
PsiDir_alpha=zeros(1,K);
% initial state Psi
PsiDir_alphasum=psi(sum(hmm.Dir_alpha,2));
% initial state KL-div
KLdiv=dirichlet_kl(hmm.Dir_alpha,hmm.prior.Dir_alpha);
for l=1:K,
    % KL-divergence for transition prob
    KLdiv=[KLdiv dirichlet_kl(hmm.Dir2d_alpha(l,:),hmm.prior.Dir2d_alpha(l,:))];
    % initial state Psi(alpha)
    PsiDir_alpha(l)=psi(hmm.Dir_alpha(l));
    for in=1:length(T);
        j = sum(T(1:in-1)) - order*(in-1) + 1;
        avLL=avLL+Gamma(j,l)*(PsiDir_alpha(l)-PsiDir_alphasum);
    end
    PsiDir2d_alphasum=psi(sum(hmm.Dir2d_alpha(l,:),2));
    for k=1:K,
        PsiDir2d_alpha(l,k)=psi(hmm.Dir2d_alpha(l,k));
        avLL=avLL+sum(Xi(:,l,k),1)*(PsiDir2d_alpha(l,k)-PsiDir2d_alphasum);
    end;
end;

XX = []; Y = [];
for in=1:N
    t0 = sum(T(1:in-1)); s0 = t0 - order*(in-1);
    XX0 = zeros(T(in)-order,length(orders)*ndim);
    for i=1:length(orders)
        o = orders(i);
        XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+order-o+1:t0+T(in)-o,:);
    end;
    XX = [XX; XX0];
    Y = [Y; residuals(s0+1:s0+T(in)-order,:)];
end

ltpi = ndim/2*log(2*pi);

if strcmp(hmm.train.covtype,'uniquediag') || strcmp(hmm.train.covtype,'uniquefull')
    if strcmp(hmm.train.covtype,'uniquediag')
        OmegaKL = 0;
        for n=1:ndim
            OmegaKL = OmegaKL + gamma_kl(hmm.Omega.Gam_shape,hmm.prior.Omega.Gam_shape, ...
                hmm.Omega.Gam_rate(n),hmm.prior.Omega.Gam_rate(n));
        end;
    else
        OmegaKL = wishart_kl(hmm.Omega.Gam_rate,hmm.prior.Omega.Gam_rate, ...
            hmm.Omega.Gam_shape,hmm.prior.Omega.Gam_shape);
    end
    KLdiv=[KLdiv OmegaKL];
end

switch hmm.train.covtype,
    case 'uniquediag'
        ldetWishB=0;
        PsiWish_alphasum=0;
        for n=1:ndim,
            ldetWishB=ldetWishB+0.5*log(hmm.Omega.Gam_rate(n));
            PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hmm.Omega.Gam_shape);
        end;
        C = hmm.Omega.Gam_shape ./ hmm.Omega.Gam_rate;
        avLL=avLL+ Tres*(-ltpi-ldetWishB+PsiWish_alphasum);
    case 'uniquefull'
        ldetWishB=0.5*logdet(hmm.Omega.Gam_rate);
        PsiWish_alphasum=0;
        for d=1:ndim,
            PsiWish_alphasum=PsiWish_alphasum+psi(hmm.Omega.Gam_shape/2+0.5-d/2);  % /2 ??
        end;
        PsiWish_alphasum=PsiWish_alphasum*0.5;
        C = hmm.Omega.Gam_shape * hmm.Omega.Gam_irate;
        avLL=avLL+ Tres*(-ltpi-ldetWishB+PsiWish_alphasum);
end;

for k=1:K,
    hs=hmm.state(k);		% for ease of referencing
    pr=hmm.state(k).prior;
    
    switch hmm.train.covtype,
        case 'diag'
            ldetWishB=0;
            PsiWish_alphasum=0;
            for n=1:ndim,
                ldetWishB=ldetWishB+0.5*log(hs.Omega.Gam_rate(n));
                PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hs.Omega.Gam_shape);
            end;
            C = hs.Omega.Gam_shape ./ hs.Omega.Gam_rate;
            avLL=avLL+ Gammasum(k)*(-ltpi-ldetWishB+PsiWish_alphasum);
        case 'full'
            ldetWishB=0.5*logdet(hs.Omega.Gam_rate);
            PsiWish_alphasum=0;
            for d=1:ndim,
                PsiWish_alphasum=PsiWish_alphasum+psi(hs.Omega.Gam_shape/2+0.5-d/2);  % /2 ??
            end;
            PsiWish_alphasum=PsiWish_alphasum*0.5;
            C = hs.Omega.Gam_shape * hs.Omega.Gam_irate;
            avLL=avLL+ Gammasum(k)*(-ltpi-ldetWishB+PsiWish_alphasum);
    end;
    
    meand =  XX * hs.W.Mu_W;
    d = residuals - meand;
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
        Cd =  repmat(C',1,Tres) .* d';
    else
        Cd = C * d';
    end
    
    dist=zeros(Tres,1);
    for n=1:ndim,
        dist=dist-0.5*d(:,n).*Cd(n,:)';
    end
    NormWishtrace=zeros(Tres,1);
    
    if order>0
        switch hmm.train.covtype,
            case 'diag'
                for n=1:ndim,
                    NormWishtrace = NormWishtrace + 0.5 * C(n) * sum( (XX * permute(hs.W.S_W(n,:,:),[2 3 1])) .* XX, 2);
                end;
                
            case 'uniquediag'
                for n=1:ndim,
                    NormWishtrace = NormWishtrace + 0.5 * C(n) * sum( (XX * permute(hs.W.S_W(n,:,:),[2 3 1])) .* XX, 2);
                end;
                
            case 'full'
                for n1=1:ndim
                    for n2=n1:ndim
                        index1 = (0:length(orders)*ndim-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim-1) * ndim + n2;
                        NormWishtrace = NormWishtrace + 0.5 * C(n1,n2) * sum( (XX * hs.W.S_W(index1,index2)) .* XX, 2);
                    end
                end
                
            case 'uniquefull'
                for n1=1:ndim
                    for n2=n1:ndim
                        index1 = (0:length(orders)*ndim-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim-1) * ndim + n2;
                        NormWishtrace = NormWishtrace + 0.5 * C(n1,n2) * sum( (XX * hs.W.S_W(index1,index2)) .* XX, 2);
                    end
                end
        end
    end
    
    avLL = avLL + sum(Gamma(:,k).*(dist - NormWishtrace));
    
    WKL = 0;
    if order>0
        if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
            alphamat = repmat( (hs.alpha.Gam_shape ./  hs.alpha.Gam_rate), ndim, 1);
            for n=1:ndim
                prior_prec = diag( repmat(hs.sigma.Gam_shape(n,:) ./  hs.sigma.Gam_rate(n,:), 1, length(orders)) .* alphamat(:)' );
                prior_var = inv(prior_prec);
                WKL = WKL + gauss_kl(hs.W.Mu_W(:,n),zeros(ndim*length(orders),1), permute(hs.W.S_W(n,:,:),[2 3 1]), prior_var);
            end;
        else
            alphaterm = repmat( (hs.alpha.Gam_shape ./  hs.alpha.Gam_rate), ndim^2, 1);
            alphaterm = alphaterm(:);
            sigmaterm = (hmm.state(k).sigma.Gam_shape ./ hmm.state(k).sigma.Gam_rate);
            sigmaterm = repmat(sigmaterm(:), length(orders), 1);
            prior_prec = diag(alphaterm .* sigmaterm);
            prior_var = inv(prior_prec);
            mu_w = hs.W.Mu_W';
            mu_w = mu_w(:);
            WKL = gauss_kl(mu_w,zeros(ndim*ndim*length(orders),1), hs.W.S_W, prior_var);
        end
    end
    
    OmegaKL = 0;
    switch hmm.train.covtype
        case 'diag'
            OmegaKL = 0;
            for n=1:ndim
                OmegaKL = OmegaKL + gamma_kl(hs.Omega.Gam_shape,hs.prior.Omega.Gam_shape, ...
                    hs.Omega.Gam_rate(n),hs.prior.Omega.Gam_rate(n));
            end;
        case 'full'
            OmegaKL = wishart_kl(hs.Omega.Gam_rate,hs.prior.Omega.Gam_rate, ...
                hs.Omega.Gam_shape,hs.prior.Omega.Gam_shape);
    end
    
    sigmaKL = 0;
    for n1=1:ndim
        for n2=n1:ndim
            sigmaKL = sigmaKL + gamma_kl(hs.sigma.Gam_shape(n1,n2),pr.sigma.Gam_shape(n1,n2), ...
                hs.sigma.Gam_rate(n1,n2),pr.sigma.Gam_rate(n1,n2));
        end;
    end;
    alphaKL = 0;
    for i=1:length(orders)
        alphaKL = alphaKL + gamma_kl(hs.alpha.Gam_shape,pr.alpha.Gam_shape, ...
            hs.alpha.Gam_rate(i),pr.alpha.Gam_rate(i));
    end;
    
    KLdiv=[KLdiv OmegaKL sigmaKL alphaKL WKL];
end;


FrEn=[Entr -avLL +KLdiv];
