function [residuals,W] = getresiduals (X,T,order,orderGL,timelag,lambdaGL,W)
%
% Compute residuals - useful when using a global model to remove global trends
%
% INPUT
% data          observations - a struct with X (time series) and C (classes)
% T             length of series
% hmm           hmm data structure
%
% OUTPUT
% residuals     
%
% Author: Diego Vidaurre, OHBA, University of Oxford

N = length(T);

if order > 0, order = 1:timelag:order; order = order(end); end
if orderGL > 0, orderGL = 1:timelag:orderGL; orderGL = orderGL(end); end

if orderGL>0
    if nargin<7 || isempty(W),
        [W,~,~,r] = mlmar(X,T,orderGL,timelag,lambdaGL);
    else
        [~,~,~,r] = mlmar(X,T,orderGL,timelag,lambdaGL,W);
    end
    residuals = [];
    for in=1:N
        t0 = sum(T(1:in-1)) - orderGL*(in-1) ;
        t1 = t0 + T(in) - orderGL;  
        residuals = [residuals; r(t0+1+order-orderGL : t1 ,:)];
    end
else
    W = [];
    residuals = [];
    for in=1:N
        t0 = sum(T(1:in-1));
        residuals = [residuals; X(t0+order+1:t0+T(in),:)];
    end
end
