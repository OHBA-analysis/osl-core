function [mcverror,cverror] = cvhmmmar (data,T,options)
%
% Obtains the cross-validated sum of prediction quadratic errors.
%
% INPUT
% data          observations, either a struct with X (time series) and C (classes, optional)
%                             or just a matrix containing the time series
% T             length of series
% options       structure with the training options - see documentation
%
% OUTPUT
% mcverror      the averaged cross-validated mean squared error
% cverror       the mean squared error per fold
%
% Author: Diego Vidaurre, OHBA, University of Oxford


if ~isfield(options,'K'), error('K was not specified'); end
if ~isfield(options,'cyc'), options.cyc = 1000; end
if ~isfield(options,'tol'), options.tol = 0.001; end
if ~isfield(options,'initcyc'), options.initcyc = 5; end
if ~isfield(options,'Gamma'), options.Gamma = []; end
if ~isfield(options,'orderGL'), options.orderGL = 0; end
if ~isfield(options,'lambdaGL'), options.lambdaGL = 0; end
if ~isfield(options,'timelag'), options.timelag = 1; end
if ~isfield(options,'whitening'), options.whitening = 0; end
if ~isfield(options,'embeddedlags'), options.embeddedlags = 0; end
if ~isfield(options,'cvfolds'), options.cvfolds = length(T); end
if ~isfield(options,'cvrep'), options.cvrep = 1; end
if ~isfield(options,'covtype'), options.covtype = 'full'; end
if ~isfield(options,'verbose'), options.verbose = 1; end

if ~isstruct(data), data = struct('X',data); end
if ~isfield(data,'C'), data.C = NaN(size(data.X,1),options.K); end
if options.K<2, error('Not enough states, please increase K'); end
if ~isfield(options,'order'), error('You need to specify a MAR order'); end
if length(options.cvfolds)>1 && length(options.cvfolds)~=length(T), error('Incorrect assigment of trials to folds'); end
if length(options.cvfolds)>1 && ~isempty(options.Gamma), error('Set options.Gamma=[] for cross-validating'); end
if length(options.cvfolds)==1 && options.cvfolds==0, error('Set options.cvfolds to a positive integer'); end
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

train = struct();
train.K = options.K;
train.order=options.order;
train.covtype = options.covtype; 
train.orderGL = options.orderGL;
train.lambdaGL = options.lambdaGL;
train.timelag = options.timelag;
train.embeddedlags = options.embeddedlags;
train.cyc = options.cyc; 
train.tol = options.tol; 
train.verbose = options.verbose; 
verbose = options.verbose;

mcverror = 0;
if length(options.cvfolds)==1,
    options.cvfolds = crossvalind('Kfold', length(T), options.cvfolds);
end
nfolds = max(options.cvfolds);
if train.order > 0, order = 1:train.timelag:train.order; order = order(end);
else order = 0; end
Ttotal = 0;
cverror = zeros(nfolds,1);
for fold=1:nfolds
    datatr.X = []; datatr.C = []; Ttr = [];
    datate.X = []; datate.C = []; Tte = []; Gammate = []; test = [];
    % build fold
    for i=1:length(T)
        t0 = sum(T(1:(i-1)))+1; t1 = sum(T(1:i)); Ti = t1-t0+1;
        if options.cvfolds(i)==fold % in testing
            datate.X = [datate.X; data.X(t0:t1,:)];
            datate.C = [datate.C; data.C(t0:t1,:)];
            Tte = [Tte Ti];
            Gammate = [Gammate; data.C(t0+order:t1,:)];
            test = [test; ones(Ti,1)];
        else % in training
            datatr.X = [datatr.X; data.X(t0:t1,:)];
            datatr.C = [datatr.C; data.C(t0:t1,:)];
            Ttr = [Ttr Ti];
        end
    end
    Ttotal = Ttotal + sum(Tte) - length(Tte)*order;
    
    if options.whitening>0
        mu = mean(datatr.X);
        datatr.X = bsxfun(@minus, datatr.X, mu);
        datate.X = bsxfun(@minus, datate.X, mu);
        [V,D] = svd(datatr.X'*datatr.X);
        A = sqrt(size(datatr.X,1)-1)*V*sqrtm(inv(D + eye(size(D))*0.00001))*V';
        datatr.X = datatr.X*A;
        datate.X = datate.X*A;
        iA = pinv(A);
    end
        
    Fe = Inf;
    for it=1:options.cvrep
        if verbose, fprintf('CV fold %d, repetition %d \n',fold,it); end
        % init Gamma
        options.Gamma = [];
        if options.initcyc>0
            options.nu = sum(Ttr)/200; options.diagcov = 0;
            options.cyc0 = 100; options.verbose = 0;
            options.Gamma = em_init(datatr,Ttr,options);
        else
            options.Gamma = [];
            order = 1:options.timelag:options.order; order = order(end);
            for in=1:length(T)
                gamma = rand(T(in)-order,options.K);
                options.Gamma = [options.Gamma; gamma ./ repmat(sum(gamma,2),1,options.K)];
            end
        end
        % train
        hmmtr=struct('train',struct());
        hmmtr.K = options.K; hmmtr.train = train; hmmtr.train.verbose = 0;
        hmmtr=hmmhsinit(hmmtr);
        [hmmtr,residualstr,W0tr]=obsinit(datatr,Ttr,hmmtr,options.Gamma);
        [hmmtr,~,~,fe,actstates] = hmmtrain(datatr,Ttr,hmmtr,options.Gamma,residualstr); fe = fe(end);
        % test
        if fe<Fe,
            Fe = fe; 
            residualste = getresiduals(datate.X,Tte,hmmtr.train.order,hmmtr.train.orderGL,hmmtr.train.timelag,hmmtr.train.lambdaGL,W0tr);
            cverror(fold) = mean(hmmerror(datate.X,Tte,hmmtr,Gammate,test,residualste,actstates));
        end
    end
    
    mcverror = mcverror + (sum(Tte) - length(Tte)*order) * cverror(fold);
end
mcverror = mcverror / Ttotal;
