function Gamma = hmmmar_init(data,T,options)
%
% Initialise the hidden Markov chain using HMM-MAR
%
% INPUT
% data      observations, a struct with X (time series) and C (classes, optional)
% T         length of observation sequence
% options,  structure with the training options  
%
% OUTPUT
% Gamma     p(state given X)
%
% Author: Diego Vidaurre, University of Oxford

[orders,order] = formorders(options.order,options.orderoffset,options.timelag,options.exptimelag);
Sind = formindexes(orders,options.S);

fehist = Inf;
for it=1:options.initrep
    options.Gamma = [];
    for in=1:length(T)
        gamma = rand(T(in)-order,options.K);
        options.Gamma = [options.Gamma; gamma ./ repmat(sum(gamma,2),1,options.K)];
    end
    hmm0=struct('train',struct());
    hmm0.K = options.K;
    hmm0.train = options; 
    hmm0.train.Sind = Sind; 
    hmm0.train.cyc = hmm0.train.initcyc;
    hmm0.train.verbose = 0;
    hmm0=hmmhsinit(hmm0);
    [hmm0,residuals0]=obsinit(data,T,hmm0,options.Gamma);
    [~,Gamma0,~,fehist0] = hmmtrain(data,T,hmm0,options.Gamma,residuals0);
    if options.verbose,
        fprintf('Init run %d, Free Energy %f \n',it,fehist0(end));
    end
    if fehist0(end)<fehist(end),
        fehist = fehist0; Gamma = Gamma0; 
    end
end

end