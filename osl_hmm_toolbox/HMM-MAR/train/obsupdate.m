function [hmm] = obsupdate (X,T,Gamma,hmm,residuals)
%
% Update observation model
%
% INPUT
% X             observations
% T             length of series
% Gamma         p(state given X)
% hmm           hmm data structure
% residuals     in case we train on residuals, the value of those.
%
% OUTPUT
% hmm           estimated HMMMAR model
%
% Author: Diego Vidaurre, OHBA, University of Oxford

[orders,order] = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);
Tres = sum(T) - length(T)*order;
ndim = size(X,2);
K=hmm.K;

obs_tol = 0.00005;
obs_maxit = 1; %20;
mean_change = Inf;
obs_it = 1;

% Some stuff that will be later used
Gammasum = sum(Gamma);
XX = formautoregr(X,T,orders,order,hmm.train.zeromean);
XXGXX = zeros(size(XX,2),size(XX,2),K);
for k=1:K
    XXGXX(:,:,k) = (XX' .* repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim*length(orders),1)) * XX;
end

while mean_change>obs_tol && obs_it<=obs_maxit,
    
    last_state = hmm.state;
        
    %%% W
    XW = zeros(K,Tres,ndim);
    if order>0 || hmm.train.zeromean==0
        [hmm,XW] = updateW(hmm,Gamma,residuals,orders,XX,XXGXX);
    end
       
    %%% Omega
    hmm = updateOmega(hmm,Gamma,Gammasum,residuals,orders,Tres,XX,XXGXX,XW);
    
    %%% sigma - channel x channel coefficients
    %%% alpha - one per order
    if order>0 
        hmm = updateSigma(hmm,orders);
        hmm = updateAlpha(hmm,orders);
    end
    
    %%% termination conditions
    obs_it = obs_it + 1;
    mean_changew = 0;
    for k=1:K
        mean_changew = mean_changew + sum(sum(abs(last_state(k).W.Mu_W - hmm.state(k).W.Mu_W))) / length(orders) / sum(hmm.train.S(:)) / K;
    end;
    mean_change = mean_changew;
end;

end
