function [e,explained,response]=hmmerror(X,T,hmm,Gamma,test,residuals,actstates)
%
% Computes the error and explained variance for the given data points  
%
% X         Observations
% T         Number of time points for each time series
% hmm       hmm data structure
% Gamma     probability of current state cond. on data - 
%           inference is run for time points with Gamma=NaN,  
% test      sum(T) x 1 vector indicating in which time points the error is to be computed;
%           common use will be to set test(t) = 1 when Gamma is NaN at this point
% residuals     in case we train on residuals, the value of those.
% actstates     Kx1 vector indicating which states were effectively used in the training
%
% e             mean quadratic error
% explained     Percentage of explained variance
% response      mean of the predictive response distribution for test(t)=1
%
% Author: Diego Vidaurre, OHBA, University of Oxford

N = length(T); 
if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end); 
else order = 0; end 

if nargin<6 || isempty(residuals),
   residuals = getresiduals(X,T,hmm.train.order,hmm.train.orderGL,hmm.train.timelag,hmm.train.lambdaGL); 
end
if nargin<7,
    actstates = ones(hmm.K,1); 
end

Y = []; 
te = [];
for in=1:N
    t0 = sum(T(1:in-1));  
    Y = [Y; X(t0+1+order:t0+T(in),:)];
    te = [te; test(t0+1+order:t0+T(in))];
end

response = hmmpred(X,T,hmm,Gamma,residuals,actstates); response = response(te==1,:);
e = (mean(residuals(te==1,:) - response).^2) ./ var(residuals(te==1,:));
explained = 1 - e ./ mean( (residuals(te==1,:) - repmat(mean(residuals(te==1,:)),sum(te),1) ).^2 );

