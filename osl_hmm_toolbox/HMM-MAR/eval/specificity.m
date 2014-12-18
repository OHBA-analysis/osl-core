function sp = specificity(X,T,hmm,Gamma)
%
% Get a measure of how specific or distinct from each other are the states 
% 
% INPUT
% X             a matrix containing the time series
% T             length of series
% hmm           the HMM-MAR structure 
% Gamma         Time courses of the states probabilities given data
%
% OUTPUT
% sp            Value of the specificity measure
%
% Author: Diego Vidaurre, OHBA, University of Oxford

order = hmm.train.order;
timelag = hmm.train.timelag;
hmm.orderGL = 0;

var0 = sum( (X - repmat(mean(X),size(X,1),1)).^2 );

% get global residuals
[~,~,~,r0] = mlmar(X,T,order,timelag);
e0 = sum(sum( r0.^2 ) ./ var0);

% get local residuals
[~,~,gcovm] = mlhmmmar(X,T,hmm,Gamma);
nsamples = sum(T) - length(T)*order;
e2 = sum( (diag(gcovm)' * nsamples) ./ var0);
sp = e2/e0;





