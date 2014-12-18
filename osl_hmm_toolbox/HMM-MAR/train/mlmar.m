function [W,covm,pred,residuals,cod] = mlmar (X,T,order,timelag,lambda,W)
%
% Estimates maximum-likelihood (ML) MAR model
% 
% INPUT
% X             observations
% T             length of series
% order         order of the ML MAR estimations
% timelag       time separation between lags
% W             MAR parameters - not to be computed if they are supplied 
% lambda        L2-regularisation parameter (ridge estimation)
%
% OUTPUT
% W             ML MAR parameters, (orderxN) by N
% covm          ML covariance matrix of the error 
% pred          predicted values
% residuals     fit residuals
% cod           coefficient of determination
%
% Author: Diego Vidaurre, OHBA, University of Oxford

ndim = size(X,2);
N = length(T);
if nargin<5, lambda=0; end
if nargin<4, timelag=1; end


XX = []; Y = [];
orders = 1:timelag:order; %orders = orders(2:end);
for in=1:N
    if in==1, t0 = 0;
    else t0 = sum(T(1:in-1));
    end
    if order>0
        XX0 = zeros(T(in)-orders(end),length(orders)*ndim);
        for i=1:length(orders)
            o = orders(i); 
            XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+orders(end)-o+1:t0+T(in)-o,:);
        end;
        XX = [XX; XX0];
    end
    Y = [Y; X(t0+orders(end)+1:t0+T(in),:)];
end


if nargin<6 || isempty(W)
    if order>0
        W = (XX' * XX + lambda * eye(size(XX,2))) \ XX' * Y;
        pred = XX * W;
    else
        pred = zeros(sum(T),ndim);
    end
else
    pred = XX * W;
end


residuals = Y - pred;
SSE = sum ( residuals.^2 );
covm = (residuals' * residuals ) / (size(Y,1)-1) ;
cod = 1 - SSE ./ sum( (Y - repmat( mean(Y),size(Y,1),1) ).^2 );


