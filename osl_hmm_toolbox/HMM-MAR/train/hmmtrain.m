function [hmm,Gamma,Xi,fehist,actstates]=hmmtrain(data,T,hmm,Gamma,residuals)
%
% Train Hidden Markov Model using using Variational Framework
%
% INPUTS:
%
% data          observations - a struct with X (time series) and C (classes)
% T             Number of time points for each time series
% hmm           hmm structure with options specified in hmm.train
% Gamma         Initial state courses
% residuals     in case we train on residuals, the value of those.
%
% OUTPUTS
% hmm           estimated HMMMAR model 
% Gamma         estimated p(state | data)
% Xi            joint probability of past and future states conditioned on data
% fehist        historic of the free energies across iterations
% knocked       states knocked out by the Bayesian inference
%
% hmm.Pi          - intial state probability
% hmm.P           - state transition matrix
% hmm.state(k).$$ - whatever parameters there are in the observation model
%
% Author: Diego Vidaurre, OHBA, University of Oxford

N = length(T);
ndim = size(data.X,2);
if hmm.train.order > 0, order = 1:hmm.train.timelag:hmm.train.order; order = order(end); 
else order = 0; end
K = hmm.train.K;
cyc = hmm.train.cyc;
tol = hmm.train.tol;

FrEn=0;
fehist=[];
actstates = ones(1,K);

for cycle=1:cyc
    
    %%%% E step
    if hmm.K>1 || cycle==1
        [Gamma,Gammasum,Xi]=hsinference(data,T,hmm,residuals);
        if size(Gammasum,1)>1, Gammasum = sum(Gammasum); end
        if (hmm.K>1 && any(round(Gammasum) == sum(T)-N*order))
            fprintf('cycle %i: All the points collapsed in one state \n',cycle)
            break
        end
    end
    as1 = find(actstates==1);
    [as,hmm,Gamma,Xi] = getactivestates(data.X,hmm,Gamma,Xi);
    actstates(as1(as==0)) = 0;
    
    %%%% Free energy computation
    oldFrEn=FrEn;
    FrEn=sum(evalfreeenergy(data.X,T,Gamma,Xi,hmm,residuals));
    fehist=[fehist; FrEn];
    if hmm.train.verbose
        fprintf('cycle %i free energy = %f, %f \n',cycle,FrEn,(FrEn - oldFrEn) /oldFrEn*100);
    end
    if cycle>2
        if ((oldFrEn - FrEn)/oldFrEn*100) < tol
            break;
        end;
    end;
    
    %%%% M STEP
    % Observation model
    hmm=obsupdate(data.X,T,Gamma,hmm,residuals);
    
    % transition matrices and initial state
    hmm=hsupdate(Xi,Gamma,T,hmm);
    
end

if hmm.train.verbose
    fprintf('Model: %d kernels, %d dimension(s), %d data samples \n',K,ndim,sum(T));
end

return;

