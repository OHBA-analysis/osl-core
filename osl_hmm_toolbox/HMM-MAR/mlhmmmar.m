function [hmm,pred,gcovm] = mlhmmmar (X,T,hmm,Gamma)
% Given the state time courses estimation, does a last estimation of each MAR using (local) ML
%
% INPUT
% X             observations
% T             length of series
% hmm           HMM-MAR structure
% Gamma         p(state given X) - has to be fully defined

%
% OUTPUT
% hmm           HMM-MAR structure with the coefficients and covariance matrices updated to 
%                   follow a maximum-likelihood estimation
% pred          predicted response
% gcovm         covariance matrix of the error for the entire model
%
% Author: Diego Vidaurre, OHBA, University of Oxford

order = hmm.train.order;
orderGL = hmm.train.orderGL;
lambdaGL = hmm.train.lambdaGL;
timelag = hmm.train.timelag;

if order > 0, orders = 1:timelag:order; order = orders(end); 
else orders = []; end 
if orderGL > 0, ordersGL = 1:timelag:orderGL; orderGL = ordersGL(end); 
else ordersGL = []; end

K = size(Gamma,2);
ndim = size(X,2);
N = length(T);
 
if orderGL>0
    [~,~,~,R] = mlmar(X,T,orderGL,timelag,lambdaGL); % Y = R;
else
    R = X;
end
Y = [];  XX = []; 
%Y = X(orderGL+1:end,:);
for in=1:N
    t0 = sum(T(1:in-1)); s0 = t0 - orderGL*(in-1);
    Y = [Y; R(s0+order-orderGL+1:s0+T(in)-orderGL,:)];
    XX0 = zeros(T(in)-order,length(orders)*ndim);
    for i=1:length(orders)
        o = orders(i);
        XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+order-o+1:t0+T(in)-o,:);
    end;
    XX = [XX; XX0];
end

pred = zeros(size(Y));
for k=1:K
    if isfield(hmm.state(k).W,'S_W'), hmm.state(k).W = rmfield(hmm.state(k).W,'S_W'); end
    hmm.state(k).W.Mu_W = pinv(XX .* repmat(sqrt(Gamma(:,k)),1,size(XX,2))) * Y;
    predk = XX * hmm.state(k).W.Mu_W;
    pred = pred + repmat(Gamma(:,k),1,ndim) .* predk;
    e = Y - predk;
    if strcmp(hmm.train.covtype,'diag')
        hmm.state(k).Omega.Gam_rate = 0.5* sum( repmat(Gamma(:,k),1,ndim) .* e.^2 );
    elseif strcmp(hmm.train.covtype,'full')
        hmm.state(k).Omega.Gam_rate =  (e' .* repmat(Gamma(:,k)',ndim,1)) * e;
        hmm.state(k).Omega.Gam_irate = inv(hmm.state(k).Omega.Gam_rate);
    end
end
ge = Y - pred;
if strcmp(hmm.train.covtype,'uniquediag')
    hmm.Omega.Gam_rate = 0.5* sum( ge.^2 );
elseif strcmp(hmm.train.covtype,'uniquefull')
    hmm.Omega.Gam_rate =  (ge' .* ge);
    hmm.Omega.Gam_irate = inv(hmm.Omega.Gam_rate);
end
gcovm = 1/size(Y,1) * (ge' * ge);

