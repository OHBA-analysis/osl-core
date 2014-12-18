function [p,LL] = emclrgr(X,Y,C,K,nu,diagcov,tol,maxit,trace)
% An EM algorithm for clusterwise regression in time series.
%
% The initial cluster values are sampled according to the following stochastic process:
% 0. t = 1
% 1. A cluster k is selected by sampling from a multinomial distribution
% 2. A permanence time r is sampled from a poisson distribution of rate nu
% 3. from t to t+r-1 we set the cluster membership to k
% 4. if t<=T goto 1
%
% Input parameters
% X: inputs, TxN matrix
% Y: outputs, TxM matrix
% C: a class vector, TX1. Only those with NaN will be clusterized 
% K: number of clusters
% nu: initialisation parameter; default T/100
% diagcov: use a diagonal covariance matrix? 
% tol: an effective zero; default 0.01
% maxit: maximum number of iterations; default 100
% trace: if TRUE, shows progress
%
% Output values
% p: TxK matrix of soft cluster assignments
% LL: Log likelihood of the model
%
% Author: Diego Vidaurre, University of Oxford


[T, M] = size(Y);
N = size(X,2);
p = zeros(T,K);
K0 = size(C,2); if isnan(K0), K0 = 0; end
unknown = isnan(C(:,1));
Tnan = sum(unknown);

if nargin<5
    nu = T/100;
end
if nargin<6
    diagcov = 1;
end
if nargin<7
    tol = 0.00001;
end
if nargin<8
    maxit = 100;
end
if nargin<9
    trace = 0;
end

% random init
t=1;
while t<=T
    r = min(poissrnd(nu),T-t+1);
    c = find(mnrnd(1,ones(1,K)/K ));
    p(t:min(t+r-1,T),c) = 1;
    t = t+r;
end
for k=1:K0 % overwrite what is fixed
    p(C==k,:) = 0; p(C==k,k) = 1; 
end

if N>0
    beta = zeros(N,M,K);
else
    mu = zeros(M,K);
end
if diagcov, omega = zeros(M,K);
else omega = zeros(M,M,K);
end
for k=1:K
    if N>0
        beta(:,:,k) = pinv(X .* repmat(sqrt(p(:,k)),1,N)) * Y; 
        %beta(:,:,k) = ( ((X' .* repmat(p(:,k)',N,1) ) * X)  + 0.001 * eye(size(X,2))  ) \ ( (X' .* repmat(p(:,k)',N,1) ) * Y );
        er = sqrt(repmat(p(:,k),1,M)) .* (Y - X * beta(:,:,k));
    else
        mu = sum(repmat(p(:,k),1,M) .* Y) / sum(p(:,k));
        er = sqrt(repmat(p(:,k),1,M)) .* (Y - repmat(mu,T,1));
    end
    if diagcov, omega(:,k) = sum(er.^2 ) / sum(p(:,k));
    else omega(:,:,k) = er' * er / sum(p(:,k));
    end
end
pi = ones(1,K)/K;

LL0 = -Inf;
LL = zeros(T,K);
for k = 1:K
    if N>0
        pred = X * beta(:,:,k);
    else
        pred = zeros(T,M);
    end
    if diagcov, LL(:,k) = logmvnpdf(Y,pred,diag(omega(:,k)))';
    else LL(:,k) = logmvnpdf(Y,pred,omega(:,:,k))';
    end
end
LL = sum(logsumexp(LL + repmat(log(pi),T,1),2));

it = 0;
% EM   
while LL-LL0 > tol*abs(LL) 
  
    if LL0>LL
        warning('error in the initialization - likelihood decreasing\n')
    end
    
    if trace
        fprintf('Iteration %d, LL: %f \n',it,LL)
    end
    LL0 = LL;
    
    % E step
    for k=1:K
        if N>0
            pred = X(unknown,:) * beta(:,:,k);
        else
            pred = zeros(Tnan,M);
        end
        if diagcov, p(unknown,k) = pi(k) * prod(normpdf( Y(unknown,:) , pred, repmat( sqrt(omega(:,k)'), length(unknown), 1)  ), 2) ;
        else p(unknown,k) = pi(k) * mvnpdf(  Y(unknown,:), pred, omega(:,:,k) );
        end   
    end
    p = p ./ repmat(sum(p,2),1,K);
    if any(isnan(p(:,1)))
        warning('There may be numerical precision issues');
        p(isnan(p(:,1)),:) = 1/K; 
    end
    
    % M step
    for k=1:K
        if N>0
            beta(:,:,k) = ( ((X' .* repmat(p(:,k)',N,1) ) * X) + 0.001 * eye(size(X,2)) ) \ ( (X' .* repmat(p(:,k)',N,1) ) * Y );
            er = sqrt(repmat(p(:,k),1,M)) .* (Y - X * beta(:,:,k));
        else
            mu = sum(repmat(p(:,k),1,M) .* Y) / sum(p(:,k));
            er = sqrt(repmat(p(:,k),1,M)) .* (Y - repmat(mu,T,1));
        end
        if diagcov, omega(:,k) = sum( er.^2 ) / sum(p(:,k));
        else omega(:,:,k) = er' * er / sum(p(:,k));
        end
    end
    pi = sum(p) / sum(p(:));
    
    % LL calculation
    LL = zeros(T,K);
    for k = 1:K
        if N>0
            pred = X * beta(:,:,k);
        else
            pred = zeros(T,M);
        end
        if diagcov, LL(:,k) = logmvnpdf(Y,pred,diag(omega(:,k)))';
        else LL(:,k) = logmvnpdf(Y,pred,omega(:,:,k))';
        end
    end
    LL = sum(logsumexp(LL + repmat(log(pi),T,1),2));

    %plot( p(:,1));
    
    it = it + 1;
    if it>maxit
        break
    end
    
    
end
