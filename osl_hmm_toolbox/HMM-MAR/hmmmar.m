function [hmm, Gamma, Xi, vpath, GammaInit, residuals, fehist] = hmmmar (data,T,options)
%
% Main function to train the HMM-MAR model, compute the Viterbi path and,
% if requested, obtain the cross-validated sum of prediction quadratic errors. 
% 
% INPUT
% data          observations, either a struct with X (time series) and C (classes, optional) 
%                             or just a matrix containing the time series
% T             length of series
% options       structure with the training options - see documentation
%
% OUTPUT
% hmm           estimated HMMMAR model 
% Gamma         Time courses of the states probabilities given data
% Xi            joint probability of past and future states conditioned on data
% vpath         most likely state path of hard assignments
% GammaInit     Time courses used after initialization.
% residuals     if the model is trained on the residuals, the value of those
% fehist        historic of the free energies across iterations
%
% Author: Diego Vidaurre, OHBA, University of Oxford

if ~isfield(options,'K'), error('K was not specified'); end
if ~isfield(options,'Fs'), options.Fs = 1; end
if ~isfield(options,'covtype'), options.covtype = 'full'; end
if ~isfield(options,'orderGL'), options.orderGL = 0; end
if ~isfield(options,'lambdaGL'), options.lambdaGL = 0; end
if ~isfield(options,'cyc'), options.cyc = 1000; end
if ~isfield(options,'tol'), options.tol = 0.001; end
if ~isfield(options,'initcyc'), options.initcyc = 5; end
if ~isfield(options,'Gamma'), options.Gamma = []; end

if ~isfield(options,'timelag'), options.timelag = 1; end
if ~isfield(options,'whitening'), options.whitening = 0; end
if ~isfield(options,'embeddedlags'), options.embeddedlags = 0; end
if ~isfield(options,'repetitions'), options.repetitions = 1; end
if ~isfield(options,'keepS_W'), options.keepS_W = 1; end
if ~isfield(options,'verbose'), options.verbose = 1; end

if ~isstruct(data), data = struct('X',data); end
if ~isfield(data,'C'), data.C = NaN(size(data.X,1),options.K); end
if options.K<2, error('Not enough states, please increase K'); end
if ~isfield(options,'order'), error('You need to specify a MAR order'); end
if options.K~=size(data.C,2), error('Matrix data.C should have K columns'); end


if length(options.embeddedlags)>1
   X = []; C = [];
   for in=1:length(T)
       [x, ind] = embedx(data.X(sum(T(1:in-1))+1:sum(T(1:in)),:),options.embeddedlags); X = [X; x ]; 
       c = data.C( sum(T(1:in-1))+1: sum(T(1:in)) , : ); c = c(ind,:); C = [C; c];
       T(in) = size(c,1);
   end
   data.X = X; data.C = C;
end

if options.whitening>0  
    mu = mean(data.X);
    data.X = bsxfun(@minus, data.X, mu);
    [V,D] = svd(data.X'*data.X);
    A = sqrt(size(data.X,1)-1)*V*sqrtm(inv(D + eye(size(D))*0.00001))*V';
    data.X = data.X*A;
    iA = pinv(A);
end
 
if isempty(options.Gamma), 
    if options.initcyc>0
        options.nu = sum(T)/200; options.diagcov = 0; options.cyc0 = 100;  
        options.Gamma = em_init(data,T,options);
    else
        options.Gamma = [];
        order = 1:options.timelag:options.order; order = order(end);  
        for in=1:length(T)
            gamma = rand(T(in)-order,options.K);
            options.Gamma = [options.Gamma; gamma ./ repmat(sum(gamma,2),1,options.K)];
        end
    end
    GammaInit = options.Gamma;
end

fehist = Inf; 
for it=1:options.repetitions
    hmm0=struct('train',struct());
    hmm0.K = options.K;
    hmm0.train.K = options.K;
    hmm0.train.Fs = options.Fs;
    hmm0.train.order=options.order;
    hmm0.train.covtype = options.covtype;
    hmm0.train.orderGL = options.orderGL;
    hmm0.train.lambdaGL = options.lambdaGL;
    hmm0.train.timelag = options.timelag;
    hmm0.train.embeddedlags = options.embeddedlags;
    hmm0.train.cyc = options.cyc;
    hmm0.train.tol = options.tol;
    hmm0.train.verbose = options.verbose;
    hmm0.train.whitening = options.whitening;
    if options.whitening, hmm0.train.A = A; hmm0.train.iA = iA;  end
    hmm0=hmmhsinit(hmm0);
    [hmm0,residuals0]=obsinit(data,T,hmm0,options.Gamma);
    [hmm0,Gamma0,Xi0,fehist0] = hmmtrain(data,T,hmm0,options.Gamma,residuals0);
    if fehist0(end)<fehist(end),
        fehist = fehist0; hmm = hmm0; 
        residuals = residuals0; Gamma = Gamma0; Xi = Xi0;
    end
end

vp = hmmdecode(data.X,T,hmm,residuals);
vpath=[];
for in=1:length(vp)
    vpath = [vpath; vp(in).q_star];
end
if ~options.keepS_W
    for i=1:hmm.K
        hmm.state(i).W.S_W = [];
    end
end