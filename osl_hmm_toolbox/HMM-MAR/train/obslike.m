function [B] = obslike (X,hmm,residuals)
%
% Evaluate likelihood of data given observation model
%
% INPUT
% X          N by ndim data matrix
% hmm        hmm data structure
% residuals  in case we train on residuals, the value of those.
%
% OUTPUT
% B          Likelihood of N data points
%
% Author: Diego Vidaurre, OHBA, University of Oxford


[T,ndim]=size(X);
K=hmm.K;
if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end);
else orders = []; order = 0; end

B=zeros(T,K);

ltpi= ndim/2 * log(2*pi);

XX = zeros(T-order,length(orders)*ndim);
for i=1:length(orders)
    o = orders(i);
    XX(:,(1:ndim) + (i-1)*ndim) = X(order-o+1:T-o,:);
end;

switch hmm.train.covtype,
    case 'uniquediag'
        ldetWishB=0;
        PsiWish_alphasum=0;
        for n=1:ndim,
            ldetWishB=ldetWishB+0.5*log(hmm.Omega.Gam_rate(n));
            PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hmm.Omega.Gam_shape);
        end;
        icov = hmm.Omega.Gam_shape ./ hmm.Omega.Gam_rate;
    case 'uniquefull'
        ldetWishB=0.5*logdet(hmm.Omega.Gam_rate);
        PsiWish_alphasum=0;
        for d=1:ndim,
            PsiWish_alphasum=PsiWish_alphasum+psi(hmm.Omega.Gam_shape/2+0.5-d/2); % /2 ??
        end;
        PsiWish_alphasum=PsiWish_alphasum*0.5;
        icov = hmm.Omega.Gam_shape * hmm.Omega.Gam_irate;
end;

for k=1:K
    hs=hmm.state(k);
    
    switch hmm.train.covtype,
        case 'diag'
            ldetWishB=0;
            PsiWish_alphasum=0;
            for n=1:ndim,
                ldetWishB=ldetWishB+0.5*log(hs.Omega.Gam_rate(n));
                PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hs.Omega.Gam_shape);
            end;
            icov = hs.Omega.Gam_shape ./ hs.Omega.Gam_rate;
        case 'full'
            ldetWishB=0.5*logdet(hs.Omega.Gam_rate);
            PsiWish_alphasum=0;
            for d=1:ndim,
                PsiWish_alphasum=PsiWish_alphasum+psi(hs.Omega.Gam_shape/2+0.5-d/2);  % /2 ??
            end;
            PsiWish_alphasum=PsiWish_alphasum*0.5;
            icov = hs.Omega.Gam_shape * hs.Omega.Gam_irate;
    end;
    
    meand = XX * hs.W.Mu_W;
    d = residuals - meand;
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
        Cd =  repmat(icov',1,T-order) .* d';
    else
        Cd = icov * d';
    end
    
    dist=zeros(T-order,1);
    for n=1:ndim,
        dist=dist-0.5*d(:,n).*Cd(n,:)';
    end
    NormWishtrace=zeros(T-order,1);
    if order>0
        switch hmm.train.covtype,
            case 'diag'
                for n=1:ndim,
                    NormWishtrace = NormWishtrace + 0.5 * (hs.Omega.Gam_shape / hs.Omega.Gam_rate(n)) * ...
                        sum( (XX * permute(hs.W.S_W(n,:,:),[2 3 1])) .* XX, 2);
                end;
                
            case 'uniquediag'
                for n=1:ndim,
                    NormWishtrace = NormWishtrace + 0.5 * (hmm.Omega.Gam_shape / hmm.Omega.Gam_rate(n)) * ...
                        sum( (XX * permute(hs.W.S_W(n,:,:),[2 3 1])) .* XX, 2);
                end;
                
            case 'full'
                for n1=1:ndim
                    for n2=n1:ndim
                        index1 = (0:length(orders)*ndim-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim-1) * ndim + n2;
                        NormWishtrace = NormWishtrace + 0.5 * icov(n1,n2) * sum( (XX * hs.W.S_W(index1,index2)) .* XX, 2);
                    end
                end
                
            case 'uniquefull'
                for n1=1:ndim
                    for n2=n1:ndim
                        index1 = (0:length(orders)*ndim-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim-1) * ndim + n2;
                        NormWishtrace = NormWishtrace + 0.5 * icov(n1,n2) * sum( (XX * hs.W.S_W(index1,index2)) .* XX, 2);
                    end
                end
        end
    end
    
    B(order+1:T,k)= + PsiWish_alphasum - ldetWishB-ltpi + dist - NormWishtrace; %-meanS
end;
B=exp(B);
