function [actstates,hmm,Gamma,Xi] = getactivestates(X,hmm,Gamma,Xi)

ndim = size(X,2);
K = hmm.K;
if hmm.train.order > 0, order = 1:hmm.train.timelag:hmm.train.order; order = order(end); 
else order = 0; end
Gammasum = sum(Gamma);

actstates = ones(1,K);
for k=1:K
    if sum(Gamma(:,k)) < size(order*ndim,2)*10, actstates(k) = 0; end
end

Gamma = Gamma(:,actstates==1);
actstates2 = actstates;

k = 1;
k0 = 1;
while k<=K
    if actstates2(k)
        k = k + 1; k0 = k0 + 1;
        continue
    end
    K = K - 1;
    hmm.state(k) = [];
    actstates2(k) = [];
    fprintf('State %d has been knocked out with %f points - there are %d left\n',k0,Gammasum(k),K)
    Gammasum(k) = [];
    k0 = k0 + 1;
end

hmm.K = K;

hmm.Dir2d_alpha = hmm.Dir2d_alpha(find(actstates==1),find(actstates==1));
hmm.Dir_alpha = hmm.Dir_alpha(find(actstates==1));
hmm.prior.Dir2d_alpha = hmm.prior.Dir2d_alpha(find(actstates==1),find(actstates==1));
hmm.prior.Dir_alpha = hmm.prior.Dir_alpha(find(actstates==1));
hmm.P = hmm.P(find(actstates==1),find(actstates==1));
hmm.Pi = hmm.Pi(find(actstates==1));

Xi = Xi(:,find(actstates),find(actstates));

end
