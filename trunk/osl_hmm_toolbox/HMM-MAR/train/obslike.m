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
[orders,order] = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);

B=zeros(T,K);
ltpi= ndim/2 * log(2*pi) - ndim/2; 
XX = formautoregr(X,size(X,1),orders,order,hmm.train.zeromean);
Sind = hmm.train.Sind==1; 
if ~hmm.train.zeromean, Sind = [true(1,size(X,2)); Sind]; end
S = hmm.train.S==1;
regressed = sum(S,1)>0;

switch hmm.train.covtype,
    case 'uniquediag'
        ldetWishB=0;
        PsiWish_alphasum=0;
        for n=1:ndim,
            if ~regressed(n), continue; end
            ldetWishB=ldetWishB+0.5*log(hmm.Omega.Gam_rate(n));
            PsiWish_alphasum=PsiWish_alphasum+0.5*psi(hmm.Omega.Gam_shape);
        end;
        C = hmm.Omega.Gam_shape ./ hmm.Omega.Gam_rate;
    case 'uniquefull'
        ldetWishB=0.5*logdet(hmm.Omega.Gam_rate(regressed,regressed));
        PsiWish_alphasum=0;
        for n=1:sum(regressed),
            PsiWish_alphasum=PsiWish_alphasum+psi(hmm.Omega.Gam_shape/2+0.5-n/2); 
        end;
        PsiWish_alphasum=PsiWish_alphasum*0.5;
        C = hmm.Omega.Gam_shape * hmm.Omega.Gam_irate;
end;

for k=1:K
    hs=hmm.state(k);
    
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
        case 'full'
            ldetWishB=0.5*logdet(hs.Omega.Gam_rate(regressed,regressed));
            PsiWish_alphasum=0;
            for n=1:sum(regressed),
                PsiWish_alphasum=PsiWish_alphasum+psi(hs.Omega.Gam_shape/2+0.5-n/2);  
            end;
            PsiWish_alphasum=PsiWish_alphasum*0.5;
            C = hs.Omega.Gam_shape * hs.Omega.Gam_irate;
    end;

    meand = zeros(size(XX,1),ndim);
    if order>0 || ~hmm.train.zeromean
        meand = XX * hs.W.Mu_W(:,regressed);
    end
    d = residuals(:,regressed) - meand;    
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
        Cd =  repmat(C(regressed)',1,T-order) .* d';
    else
        Cd = C(regressed,regressed) * d';
    end
    
    dist=zeros(T-order,1);
    for n=1:sum(regressed),
        dist=dist-0.5*d(:,n).*Cd(n,:)';
    end
    
    NormWishtrace=zeros(T-order,1);
    if order>0 || ~hmm.train.zeromean
        switch hmm.train.covtype,
            case {'diag','uniquediag'}
                for n=1:ndim,
                    if ~regressed(n), continue; end
                    NormWishtrace = NormWishtrace + 0.5 * C(n) * ...
                        sum( (XX(:,Sind(:,n)) * permute(hs.W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1])) .* XX(:,Sind(:,n)), 2);
                end;
                
            case {'full','uniquefull'}
                for n1=1:ndim
                    for n2=1:ndim
                        if ~regressed(n1) || ~regressed(n2), continue; end       
                        index1 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n1;
                        index2 = (0:length(orders)*ndim+(~hmm.train.zeromean)-1) * ndim + n2;                      
                        index1 = index1(Sind(:,n1)); index2 = index2(Sind(:,n2));
                        NormWishtrace = NormWishtrace + 0.5 * C(n1,n2) * sum( (XX * hs.W.S_W(index1,index2)) .* XX, 2);
                    end
                end
        end
    end
    
    B(order+1:T,k)= + PsiWish_alphasum - ldetWishB-ltpi + dist - NormWishtrace; %-meanS
end;
B=exp(B);

