function [Gamma,LL] = em_init(data,T,options)
%
% Initialise the hidden Markov chain using an EM algorithm for clusterwise regression in time series.
%
% INPUT
% data      observations, a struct with X (time series) and C (classes, optional)
% T         length of observation sequence
% options,  structure with the training options - different from HMMMAR are
%   nu            initialisation parameter; default T/200
%   diagcov   use a diagonal covariance matrix?
%   cyc0     maximum number of iterations; default 100
%
% OUTPUT
% Gamma     p(state given X)
% LL        the final model log-likelihood
%
% Author: Diego Vidaurre, University of Oxford

if options.order > 0, orders = 1:options.timelag:options.order; order = orders(end); 
else orders = []; order = 0; end 
if options.orderGL > 0, orderGL = 1:options.timelag:options.orderGL; orderGL = orderGL(end);
else orderGL = 0; end 
N = length(T);
ndim = size(data.X,2);
Gamma = [];

Y = []; C = [];
if orderGL>0
    [~,~,~,r] = mlmar(data.X,T,orderGL,options.timelag,options.lambdaGL);
    for in=1:N
        t0 = sum(T(1:in-1)) - orderGL*(in-1) ;
        t1 = t0 + T(in) - orderGL;  
        Y = [Y; r(t0+1+order-orderGL:t1,:)];
        C = [C; data.C(t0+1+order:t0+T(in),:)];
    end
else
    for in=1:N
        t0 = sum(T(1:in-1)); 
        Y = [Y; data.X(t0+1+order:t0+T(in),:)];
        C = [C; data.C(t0+1+order:t0+T(in),:)];
    end
end
    
XX = []; 
for in=1:N
    t0 = sum(T(1:in-1));
    XX0 = zeros(T(in)-order,length(orders)*ndim);
    for i=1:length(orders)
        o = orders(i);
        XX0(:,(1:ndim) + (i-1)*ndim) = data.X(t0+orders(end)-o+1:t0+T(in)-o,:);
    end;
    XX = [XX; XX0];
end

LL = -Inf; 
for n=1:options.initcyc
    while 1
        [gamma,ll] = emclrgr(XX,Y,C,options.K,options.nu,options.diagcov,options.tol,options.cyc0,0);
        %if order>0, [gamma,ll] = emclrgr(XX,Y,K,nu,diagcov,tol,maxit,0);
        %else [gamma,ll] = emmixg(X,K,nu,diagcov,tol,maxit,0);
        %end
        if (~isinf(ll)), break;
        else fprintf('Overflow \n');
        end
    end
    
    %er = zeros(size(Y));
    %for k=1:size(gamma,2)
    %    beta = pinv(XX .* repmat(sqrt(gamma(:,k)),1,size(XX,2) )) * Y;
    %    er = er + repmat(gamma(:,k),1,size(Y,2)) .* (XX*beta);
    %end
    %er = sqrt(sum(sum((er - Y).^2)));
    
    if options.verbose
        fprintf('Init run %d - %d, LL %f \n',in,n,ll);
    end
    if ll>LL
        LL = ll;
        Gamma = gamma;
        s = n;
    end
end
if options.verbose
    fprintf('%i-th was the best iteration with LL=%f \n',s,LL)
end


