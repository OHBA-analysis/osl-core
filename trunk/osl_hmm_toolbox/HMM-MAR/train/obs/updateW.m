function [hmm,XW] = updateW(hmm,Gamma,residuals,orders,XX,XXGXX)

K = length(hmm.state); ndim = size(hmm.state(1).W.Mu_W,2);
Sind = hmm.train.Sind==1; S = hmm.train.S==1;
if ~hmm.train.zeromean, Sind = [true(1,size(hmm.state(1).W.Mu_W,2)); Sind]; end
if isempty(orders), 
    order = 0; 
else
    order = orders(end);
end

if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
    for k=1:K,
        for n=1:ndim
            ndim_n = sum(S(:,n));
            if ndim_n==0 && hmm.train.zeromean==1, continue; end
            if strcmp(hmm.train.covtype,'diag')
                omega = hmm.state(k).Omega;
            else
                omega = hmm.Omega;
            end
            regterm = [];
            if hmm.train.zeromean==0
                regterm = hmm.state(k).prior.Mean.iS(n);
            end
            if order>0
                alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate), ndim_n, 1);
                regterm = [regterm; repmat(hmm.state(k).sigma.Gam_shape(S(:,n),n) ./ ...
                    hmm.state(k).sigma.Gam_rate(S(:,n),n), length(orders), 1).*alphaterm(:) ];
            end
            regterm = diag(regterm);
            prec =  regterm + (omega.Gam_shape / omega.Gam_rate(n)) * XXGXX(Sind(:,n),Sind(:,n),k);
            hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)) = inv(prec);
            hmm.state(k).W.Mu_W(Sind(:,n),n) = (( permute(hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)),[2 3 1]) * (omega.Gam_shape / omega.Gam_rate(n)) * XX(:,Sind(:,n))') .* ...
                repmat(Gamma(:,k)',sum(Sind(:,n)),1)) * residuals(:,n);
        end;
    end;
else % this only works if all(S(:)==1);  any(S(:)==0) is just not yet implemented
    for k=1:K
        mlW = (( XXGXX(:,:,k) \ XX') .* repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim*length(orders),1) * residuals)';
        regterm = [];
        if hmm.train.zeromean==0
            regterm = hmm.state(k).prior.Mean.iS;
        end
        if order>0
            sigmaterm = (hmm.state(k).sigma.Gam_shape(S) ./ hmm.state(k).sigma.Gam_rate(S));
            sigmaterm = repmat(sigmaterm, length(orders), 1);
            alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate), sum(S(:)), 1);
            alphaterm = alphaterm(:);
            regterm = [regterm; (alphaterm .* sigmaterm)];
        end
        regterm = diag(regterm);
        if strcmp(hmm.train.covtype,'full')
            omega = hmm.state(k).Omega;
        else
            omega = hmm.Omega;
        end
        prec = omega.Gam_shape * omega.Gam_irate;
        gram = kron(XXGXX(:,:,k), prec);
        hmm.state(k).W.S_W = inv(regterm + gram); 
        muW = hmm.state(k).W.S_W * gram * mlW(:);
        hmm.state(k).W.Mu_W = reshape(muW,ndim,~hmm.train.zeromean+ndim*length(orders))';
    end
end

XW = zeros(K,size(XX,1),ndim);
for k=1:K
    XW(k,:,:) = XX * hmm.state(k).W.Mu_W;
end;

end