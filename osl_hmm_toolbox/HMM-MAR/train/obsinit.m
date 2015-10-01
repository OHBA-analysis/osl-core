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
[residuals,W0] =  getresiduals(data.X,T,hmm.train.Sind,hmm.train.maxorder,hmm.train.order,...
    hmm.train.orderoffset,hmm.train.timelag,hmm.train.exptimelag,hmm.train.zeromean);
hmm=initpriors(data.X,T,hmm,residuals);
hmm=initpost(data.X,T,hmm,residuals,Gamma);

end


function [hmm] = initpriors(X,T,hmm,residuals)
% define priors

ndim = size(X,2);
rangresiduals2 = (range(residuals)/2).^2;

for k=1:hmm.K,
    if isfield(hmm.train,'state') && isfield(hmm.train.state(k),'train') && ~isempty(hmm.train.state(k).train)
        train = hmm.train.state(k).train;
    else
        train = hmm.train;
    end
    orders = formorders(train.order,train.orderoffset,train.timelag,train.exptimelag);
    if strcmp(train.covtype,'diag') || strcmp(train.covtype,'full')
        defstateprior(k)=struct('sigma',[],'alpha',[],'Omega',[],'Mean',[]);
    else
        defstateprior(k)=struct('sigma',[],'alpha',[],'Mean',[]);
    end
    if ~train.uniqueAR && isempty(train.prior)
        defstateprior(k).sigma = struct('Gam_shape',[],'Gam_rate',[]);
        defstateprior(k).sigma.Gam_shape = 0.1*ones(ndim,ndim); %+ 0.05*eye(ndim);
        defstateprior(k).sigma.Gam_rate = 0.1*ones(ndim,ndim);%  + 0.05*eye(ndim);
    end
    if ~isempty(orders) && isempty(train.prior)
        defstateprior(k).alpha = struct('Gam_shape',[],'Gam_rate',[]);
        defstateprior(k).alpha.Gam_shape = 0.1;
        defstateprior(k).alpha.Gam_rate = 0.1*ones(1,length(orders));
    end
    if ~train.zeromean, 
        defstateprior(k).Mean = struct('Mu',[],'iS',[]);
        defstateprior(k).Mean.Mu = zeros(ndim,1);
        defstateprior(k).Mean.S = rangresiduals2';
        defstateprior(k).Mean.iS = 1./rangresiduals2';
    end
    priorcov_rate = rangeerror(X,T,hmm.train.maxorder,orders,residuals);
    if strcmp(train.covtype,'full')
        defstateprior(k).Omega.Gam_rate = diag(priorcov_rate);
        defstateprior(k).Omega.Gam_shape = ndim+0.1-1;
    elseif strcmp(train.covtype,'diag')
        defstateprior(k).Omega.Gam_rate = 0.5 * priorcov_rate;
        defstateprior(k).Omega.Gam_shape = 0.5 * (ndim+0.1-1);
    end
    
end;

if strcmp(hmm.train.covtype,'uniquefull')
    hmm.prior.Omega.Gam_shape = ndim+0.1-1;
    hmm.prior.Omega.Gam_rate = diag(priorcov_rate);
elseif strcmp(hmm.train.covtype,'uniquediag')
    hmm.prior.Omega.Gam_shape = 0.5 * (ndim+0.1-1);
    hmm.prior.Omega.Gam_rate = 0.5 * priorcov_rate;
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
end

% moving the state options for convenience
for k=1:hmm.K,
    if isfield(hmm.train,'state') && isfield(hmm.train.state(k),'train') && ~isempty(hmm.train.state(k).train)
        hmm.state(k).train = hmm.train.state(k).train;
    end
end
if isfield(hmm.train,'state')
    hmm.train = rmfield(hmm.train,'state');
end
end


function [hmm] = initpost(X,T,hmm,residuals,Gamma)
% Initialising the posteriors

Tres = sum(T) - length(T)*hmm.train.maxorder;
ndim = size(X,2);
K = hmm.K;
S = hmm.train.S==1; regressed = sum(S,1)>0;

% initial random values for the states - multinomial
Gammasum = sum(Gamma);
XXGXX = cell(K,1);
setxx; % build XX and get orders
 
% W
for k=1:K
    setstateoptions; 
    hmm.state(k).W = struct('Mu_W',[],'S_W',[]);
    if order>0 || ~train.zeromean
        if train.uniqueAR || ndim==1 % it is assumed that order>0 and cov matrix is diagonal
            XY = zeros(length(orders),1);
            XGX = zeros(length(orders));
            for n=1:ndim
                ind = n:ndim:size(XX{kk},2);
                XGX = XGX + XXGXX{k}(ind,ind);
                XY = XY + (XX{kk}(:,ind)' .* repmat(Gamma(:,k)',length(ind),1)) * residuals(:,n);
            end
            if ~isempty(train.prior)
                hmm.state(k).W.S_W = inv(train.prior.iS + XGX);
                hmm.state(k).W.Mu_W = hmm.state(k).W.S_W * (XY + train.prior.iSMu); % order by 1
            else
                %hmm.state(k).W.S_W = inv(0.1 * mean(trace(XGX)) * eye(length(orders)) + XGX);
                hmm.state(k).W.S_W = inv(0.01 * eye(length(orders)) + XGX);
                hmm.state(k).W.Mu_W = hmm.state(k).W.S_W * XY; % order by 1
            end
            
        elseif strcmp(train.covtype,'uniquediag') || strcmp(train.covtype,'diag')
            hmm.state(k).W.Mu_W = zeros((~train.zeromean)+ndim*length(orders),ndim);
            hmm.state(k).W.S_W = zeros(ndim,(~train.zeromean)+ndim*length(orders),(~train.zeromean)+ndim*length(orders));
            for n=1:ndim
                ndim_n = sum(S(:,n));
                if ndim_n==0, continue; end
                hmm.state(k).W.S_W(n,Sind(:,n),Sind(:,n)) = inv( ((XX{kk}(:,Sind(:,n))' .* repmat(Gamma(:,k)',(~train.zeromean)+ndim_n*length(orders),1) ) * ...
                    XX{kk}(:,Sind(:,n))) + 0.01*eye((~train.zeromean) + ndim_n*length(orders))  )   ;
                hmm.state(k).W.Mu_W(Sind(:,n),n) = (( permute(hmm.state(k).W.S_W(1,Sind(:,n),Sind(:,n)),[2 3 1]) * XX{kk}(:,Sind(:,n))') .* ...
                    repmat(Gamma(:,k)',(~train.zeromean)+ndim_n*length(orders),1)) * residuals(:,n);
            end;
        else
            hmm.state(k).W.S_W = inv( kron(XXGXX{k},eye(ndim)) );
            hmm.state(k).W.Mu_W = (( XXGXX{k} \ XX{kk}' ) .* repmat(Gamma(:,k)',(~train.zeromean)+ndim*length(orders),1)) * residuals;
        end
    end
end;

% Omega
if strcmp(hmm.train.covtype,'uniquediag') && hmm.train.uniqueAR
    hmm.Omega.Gam_rate = hmm.prior.Omega.Gam_rate;
    for k=1:K
        if hmm.train.multipleConf, kk = k;
        else kk = 1;
        end
        XW = zeros(size(XX,1),ndim);
        for n=1:ndim
            ind = n:ndim:size(XX{kk},2);
            XW(:,n) = XX{kk}(:,ind) * hmm.state(k).W.Mu_W;
        end
        e = (residuals - XW).^2;
        hmm.Omega.Gam_rate = hmm.Omega.Gam_rate + ...
            0.5 * sum( repmat(Gamma(:,k),1,ndim) .* e );
    end;
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;    
    
elseif strcmp(hmm.train.covtype,'uniquediag')
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;
    hmm.Omega.Gam_rate = zeros(1,ndim);
    hmm.Omega.Gam_rate(regressed) = hmm.prior.Omega.Gam_rate(regressed);
    for k=1:hmm.K
        if ~isempty(hmm.state(k).W.Mu_W(:,regressed))
            e = residuals(:,regressed) - XX{kk} * hmm.state(k).W.Mu_W(:,regressed);
        else
            e = residuals(:,regressed);
        end
        hmm.Omega.Gam_rate(regressed) = hmm.Omega.Gam_rate(regressed) +  ...
            0.5 * sum( repmat(Gamma(:,k),1,sum(regressed)) .* e.^2 );
    end
    
elseif strcmp(hmm.train.covtype,'uniquefull')
    hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres;
    hmm.Omega.Gam_rate = zeros(ndim,ndim); hmm.Omega.Gam_irate = zeros(ndim,ndim);
    hmm.Omega.Gam_rate(regressed,regressed) = hmm.prior.Omega.Gam_rate(regressed,regressed);
    for k=1:hmm.K
        if ~isempty(hmm.state(k).W.Mu_W(:,regressed))
            e = residuals(:,regressed) - XX{kk} * hmm.state(k).W.Mu_W(:,regressed);
        else
            e = residuals(:,regressed);
        end
        hmm.Omega.Gam_rate(regressed,regressed) = hmm.Omega.Gam_rate(regressed,regressed) +  ...
            (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
    end
    hmm.Omega.Gam_irate(regressed,regressed) = inv(hmm.Omega.Gam_rate(regressed,regressed));

    
else % state dependent
    for k=1:K
        setstateoptions;
        if train.uniqueAR
            XW = zeros(size(XX{kk},1),ndim);
            for n=1:ndim
                ind = n:ndim:size(XX{kk},2);
                XW(:,n) = XX{kk}(:,ind) * hmm.state(k).W.Mu_W;
            end
            e = (residuals - XW).^2;
            hmm.state(k).Omega.Gam_rate = hmm.state(k).prior.Omega.Gam_rate + ...
                0.5* sum( repmat(Gamma(:,k),1,ndim) .* e );
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Gammasum(k) / 2;
                
        elseif strcmp(train.covtype,'diag')
            if ~isempty(hmm.state(k).W.Mu_W)  
                e = (residuals(:,regressed) - XX{kk} * hmm.state(k).W.Mu_W(:,regressed)).^2;
            else
                e = residuals(:,regressed).^2;
            end
            hmm.state(k).Omega.Gam_rate(regressed) = hmm.state(k).prior.Omega.Gam_rate(regressed) + ...
                sum( repmat(Gamma(:,k),1,sum(regressed)) .* e ) / 2;
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Gammasum(k) / 2;
            
        else % full
            if ~isempty(hmm.state(k).W.Mu_W)  
                e = residuals(:,regressed) - XX{kk} * hmm.state(k).W.Mu_W(:,regressed);
            else
                e = residuals(:,regressed);
            end
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + Gammasum(k);
            hmm.state(k).Omega.Gam_rate = zeros(ndim,ndim); hmm.state(k).Omega.Gam_irate = zeros(ndim,ndim);
            hmm.state(k).Omega.Gam_rate(regressed,regressed) =  hmm.state(k).prior.Omega.Gam_rate(regressed,regressed) +  ...
                (e' .* repmat(Gamma(:,k)',sum(regressed),1)) * e;
            hmm.state(k).Omega.Gam_irate(regressed,regressed) = inv(hmm.state(k).Omega.Gam_rate(regressed,regressed));
        end
    end
    
end

for k=1:K,
    if isfield(hmm.state(k),'train') && ~isempty(hmm.state(k).train), train = hmm.state(k).train;
    else train = hmm.train;
    end
    if train.order>0 && isempty(train.prior)
        hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape;
        hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    end
end

%%% sigma - channel x channel coefficients
hmm = updateSigma(hmm);

%%% alpha - one per order
hmm = updateAlpha(hmm);

end

