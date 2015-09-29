function [residuals,W] = getresiduals (X,T,Sind,maxorder,order,orderoffset,timelag,exptimelag,zeromean,W)
%
% Compute residuals - useful when using a global model to remove global trends
%
% INPUT
% X             time series
% T             length of series
%
%
% OUTPUT
% residuals
% W
%
% Author: Diego Vidaurre, OHBA, University of Oxford

if nargin<10, W = []; end

N = length(T);

[~,order] = formorders(order,orderoffset,timelag,exptimelag);

if any(Sind(:)==0) 
    if isempty(W),
        [W,~,~,residuals] = mlmar(X,T,Sind==0,maxorder,order,orderoffset,timelag,exptimelag,zeromean);
    else
        [~,~,~,residuals] = mlmar(X,T,Sind==0,maxorder,order,orderoffset,timelag,exptimelag,zeromean,W);
    end
else
    W = [];
    residuals = [];
    for in=1:N
        t0 = sum(T(1:in-1));
        residuals = [residuals; X(t0+maxorder+1:t0+T(in),:)];
    end
end
