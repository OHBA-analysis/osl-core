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

[options,data] = checkoptions(options,data,T,0);

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

[orders,order] = formorders(options.order,options.orderoffset,options.timelag,options.exptimelag);
Sind = formindexes(orders,options.S);


if isempty(options.Gamma),
    if options.K > 1
        if options.initrep>0 && strcmp(options.inittype,'HMM-MAR')
            options.Gamma = hmmmar_init(data,T,options);
        elseif options.initrep>0 && strcmp(options.inittype,'EM')
            options.nu = sum(T)/200;
            options.Gamma = em_init(data,T,options,Sind,options.zeromean);
        else
            options.Gamma = [];
            for in=1:length(T)
                gamma = rand(T(in)-order,options.K);
                options.Gamma = [options.Gamma; gamma ./ repmat(sum(gamma,2),1,options.K)];
            end
        end
    else
       options.Gamma = ones(sum(T)-length(T)*order,1); 
    end
end   

GammaInit = options.Gamma;

fehist = Inf;
for it=1:options.repetitions
    hmm0=struct('train',struct());
    hmm0.K = options.K;
    hmm0.train = options; 
    hmm0.train.Sind = Sind; 
    if options.whitening, hmm0.train.A = A; hmm0.train.iA = iA;  end
    hmm0=hmmhsinit(hmm0);
    [hmm0,residuals0]=obsinit(data,T,hmm0,options.Gamma);
    [hmm0,Gamma0,Xi0,fehist0] = hmmtrain(data,T,hmm0,options.Gamma,residuals0);
    if options.updateGamma==1 && fehist0(end)<fehist(end),
        fehist = fehist0; hmm = hmm0; 
        residuals = residuals0; Gamma = Gamma0; Xi = Xi0;
    elseif options.updateGamma==0,
        fehist = []; hmm = hmm0; 
        residuals = []; Gamma = options.Gamma; Xi = [];
    end
end

if options.updateGamma
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
else
    vpath = ones(size(Gamma,1),1);
end
hmm.train = rmfield(hmm.train,'Sind');

end