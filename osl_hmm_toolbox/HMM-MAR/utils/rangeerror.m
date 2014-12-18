function [r] = rangeerror(X,T,orders,residuals,lambda)
% estimates the range of the error

if isempty(orders),
    r = range(X);
    
else
    if nargin<5
        lambda=0.01;
    end
    ndim=size(X,2);
    XX = []; Y = [];
    for in=1:length(T)
        t0 = sum(T(1:in-1)); s0 = t0 - orders(end)*(in-1);
        XX0 = zeros(T(in)-orders(end),length(orders)*ndim);
        for i=1:length(orders)
            o = orders(i);
            XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+orders(end)-o+1:t0+T(in)-o,:);
        end;
        XX = [XX; XX0];
        Y = [Y; residuals(s0+1:s0+T(in)-orders(end),:)];
    end
    
    W = (XX' * XX + lambda*eye(length(orders)*ndim)) \ XX' * Y;
    r = Y - XX * W;
    %r = 0.5*sum(r.^2 );
    r = range(r);
    
end
