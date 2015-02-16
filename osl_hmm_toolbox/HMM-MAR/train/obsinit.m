function [hmm, residuals, W0] = obsinit (data,T,hmm,Gamma)
%
% Initialise observation model in HMM
%
% INPUT
% data          observations - a struct with X (time series) and C (classes)
% T             length of series
% Gamma         p(state given X)
% hmm           hmm data structure
%
% OUTPUT
% hmm           estimated HMMMAR model
% residuals     in case we train on residuals, the value of those
% W0            global MAR estimation
%
% Author: Diego Vidaurre, OHBA, University of Oxford

if nargin<4, Gamma = []; end
[residuals,W0] =  getresiduals(data.X,T,hmm.train.Sind,hmm.train.order,...
    hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag,hmm.train.zeromean);
hmm=initpriors(data.X,T,hmm,residuals);
hmm=initpost(data.X,T,hmm,residuals,Gamma);

end


function [hmm] = initpriors(X,T,hmm,residuals)
% define priors

orders = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);

ndim = size(X,2);
rangresiduals2 = (range(residuals)/2).^2;

for k=1:hmm.K,
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'full')
        defstateprior(k)=struct('sigma',[],'alpha',[],'Omega',[],'Mean',[]);
    else
        defstateprior(k)=struct('sigma',[],'alpha',[],'Mean',[]);
    end
    defstateprior(k).sigma = struct('Gam_shape',[],'Gam_rate',[]);
    defstateprior(k).sigma.Gam_shape = 0.1*ones(ndim,ndim); %+ 0.05*eye(ndim);
    defstateprior(k).sigma.Gam_rate = 0.1*ones(ndim,ndim);%  + 0.05*eye(ndim);
    defstateprior(k).alpha = struct('Gam_shape',[],'Gam_rate',[]);
    defstateprior(k).alpha.Gam_shape = 0.1;
    defstateprior(k).alpha.Gam_rate = 0.1*ones(1,length(orders));
    if ~hmm.train.zeromean, 
        defstateprior(k).Mean = struct('Mu',[],'iS',[]);
        defstateprior(k).Mean.Mu = zeros(ndim,1);
        defstateprior(k).Mean.S = rangresiduals2';
        defstateprior(k).Mean.iS = 1./rangresiduals2';
    end
    priorcov_rate = rangeerror(X,T,orders,residuals);
    
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'full')
        defstateprior(k).Omega.Gam_shape = ndim+0.1-1;
        if strcmp(hmm.train.covtype,'full')
            defstateprior(k).Omega.Gam_rate = diag(priorcov_rate);
        else
            defstateprior(k).Omega.Gam_rate = priorcov_rate;
        end
    end
    
end;

if strcmp(hmm.train.covtype,'uniquediag') || strcmp(hmm.train.covtype,'uniquefull')
    hmm.prior.Omega.Gam_shape = ndim+0.1-1;
    if strcmp(hmm.train.covtype,'uniquefull')
        hmm.prior.Omega.Gam_rate = diag(priorcov_rate);
    else
        hmm.prior.Omega.Gam_rate = priorcov_rate;
    end
end

% assigning default priors for observation models
if ~isfield(hmm,'state') || ~isfield(hmm.state,'prior'),
    for k=1:hmm.K,
        hmm.state(k).prior=defstateprior(k);
    end;
else
    for k=1:hmm.K,
        % prior not specified are set to default
        statepriorlist=fieldnames(defstateprior(k));
        fldname=fieldnames(hmm.state(k).prior);
        misfldname=find(~ismember(statepriorlist,fldname));
        for i=1:length(misfldname),
            priorval=getfield(defstateprior(k),statepriorlist{i});
            hmm.state(k).prior=setfield(hmm.state,k,'prior',statepriorlist{i}, ...
                priorval);
        end;
    end;
end;

end


function [hmm] = initpost(X,T,hmm,residuals,Gamma)
% Initialising the posteriors

[orders,order] = formorders(hmm.train.order,hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag);

Tres = sum(T) - length(T)*order;
ndim = size(X,2);
K = hmm.K;
Sind = hmm.train.Sind==1; S = hmm.train.S==1;
if ~hmm.train.zeromean, Sind = [true(1,size(X,2)); Sind]; end
regressed = sum(S,1)>0;

% initial random values for the states - multinomial
gammasum = sum(Gamma);

XX = formautoregr(X,T,orders,order,hmm.train.zeromean);

% W
for k=1:K,
    hmm.state(k).W = struct('Mu_W',[],'S_W',[]);
    if order>0 || ~hmm.train.zeromean
        if strcmp(hmm.train.covtype,'uniquediag') || strcmp(hmm.train.covtype,'diag')
            hmm.state(k).W.Mu_W = zeros((~hmm.train.zeromean)+ndim*length(orders),ndim);
            hmm.state(k).W.S_W = zeros(ndim,(~hmm.train.zeromean)+ndim*length(orders),(~hmm.train.zeromean)+ndim*length(orders));
            for n=1:ndim
                ndim_n = sum(S(:,n));
                if ndim_n==0, continue; end
                hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)) = inv( ((XX(:,Sind(:,n))' .* repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim_n*length(orders),1) ) * ...
                    XX(:,Sind(:,n))) + 0.01*eye((~hmm.train.zeromean) + ndim_n*length(orders))  )   ;
                hmm.state(k).W.Mu_W(Sind(:,n),n) = (( permute(hmm.state(k).W.S_W(1,Sind(:,n),Sind(:,n)),[2 3 1]) * XX(:,Sind(:,n))') .* ...
                    repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim_n*length(orders),1)) * residuals(:,n);
            end;
        else
            XXGXX = (XX' .* repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim*length(orders),1)) * XX;
            hmm.state(k).W.S_W = inv( kron(XXGXX,eye(ndim)) );
            hmm.state(k).W.Mu_W = (( XXGXX \ XX' ) .* repmat(Gamma(:,k)',(~hmm.train.zeromean)+ndim*length(orders),1)) * residuals;
        end
    end
end;

% Omega
switch hmm.train.covtype
    case 'uniquediag'
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;
        hmm.Omega.Gam_rate = zeros(1,ndim);
        hmm.Omega.Gam_rate(regressed) = hmm.prior.Omega.Gam_rate(regressed);
        for k=1:hmm.K
            if order>0
                e = residuals(:,regressed) - XX * hmm.state(k).W.Mu_W(:,regressed);
            else
                e = residuals(:,regressed);
            end
            hmm.Omega.Gam_rate(regressed) = hmm.Omega.Gam_rate(regressed) +  ...
                0.5 * sum( repmat(Gamma(:,k),1,sum(regressed)) .* e.^2 );
        end
    case 'uniquefull'
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres;
        hmm.Omega.Gam_rate = zeros(ndim,ndim); hmm.Omega.Gam_irate = zeros(ndim,ndim); 
        hmm.Omega.Gam_rate(regressed,regressed) = hmm.prior.Omega.Gam_rate(regressed,regressed);
        for k=1:hmm.K
            if order>0
                e = residuals(:,regressed) - XX * hmm.state(k).W.Mu_W(:,regressed);
            else
                e = residuals(:,regressed);
            end            
            hmm.Omega.Gam_rate(regressed,regressed) = hmm.Omega.Gam_rate(regressed,regressed) +  ...
                (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
        end
        hmm.Omega.Gam_irate(regressed,regressed) = inv(hmm.Omega.Gam_rate(regressed,regressed));
    case 'diag'
        for k=1:hmm.K
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k) / 2;
            if order>0
                e = residuals(:,regressed) - XX * hmm.state(k).W.Mu_W(:,regressed);
            else
                e = residuals(:,regressed);
            end            
            hmm.state(k).Omega.Gam_rate = zeros(1,ndim);
            hmm.state(k).Omega.Gam_rate(regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed) + ...
                0.5* sum( repmat(Gamma(:,k),1,sum(regressed)) .* e.^2 );
        end
    case 'full'
        for k=1:hmm.K
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k);
            if order>0
                e = residuals(:,regressed) - XX * hmm.state(k).W.Mu_W(:,regressed);
            else
                e = residuals(:,regressed);
            end            
            hmm.Omega.Gam_rate = zeros(ndim,ndim); hmm.state(k).Omega.Gam_irate = zeros(ndim,ndim); 
            hmm.state(k).Omega.Gam_rate(regressed,regressed) =  hmm.state(k).prior.Omega.Gam_rate(regressed,regressed) +  ...
                (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
            hmm.state(k).Omega.Gam_irate(regressed,regressed) = inv(hmm.state(k).Omega.Gam_rate(regressed,regressed));
        end
end


if order>0 
    for k=1:K,
        hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape;
        hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    end    
end
% sigma and alpha
if order>0 
    hmm = updateSigma(hmm,orders);
    hmm = updateAlpha(hmm,orders);
end


end

