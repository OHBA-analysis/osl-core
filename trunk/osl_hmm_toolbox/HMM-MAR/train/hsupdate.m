function [hmm]=hsupdate(Xi,Gamma,T,hmm)
%
% updates hidden state parameters of an HMM
%
% INPUT:
%
% Xi     probability of past and future state cond. on data
% Gamma  probability of current state cond. on data
% T      length of observation sequences
% hmm    single hmm data structure
%
% OUTPUT
% hmm    single hmm data structure with updated state model probs.
%
% Author: Diego Vidaurre, OHBA, University of Oxford


N = length(T);
K=hmm.K;
if hmm.train.order > 0, order = 1:hmm.train.timelag:hmm.train.order; order = order(end); 
else order = 0; end 
% transition matrix
sxi=squeeze(sum(Xi,1));   % counts over time

hmm.Dir2d_alpha=sxi+hmm.prior.Dir2d_alpha;
PsiSum=psi(sum(hmm.Dir2d_alpha(:),1));
for j=1:K,
    for i=1:K,
        P(j,i)=exp(psi(hmm.Dir2d_alpha(j,i))-PsiSum);
    end;
    P(j,:)=P(j,:)./sum(P(j,:));
end;
hmm.P=P;

% initial state
hmm.Dir_alpha=hmm.prior.Dir_alpha;
for in=1:N
    t = sum(T(1:in-1)) - order*(in-1) + 1;
    hmm.Dir_alpha=hmm.Dir_alpha+Gamma(t,:);
end
PsiSum=psi(sum(hmm.Dir_alpha,2));
for i=1:K,
    Pi(i)=exp(psi(hmm.Dir_alpha(i))-PsiSum);
end
hmm.Pi=hmm.Pi./sum(hmm.Pi);

