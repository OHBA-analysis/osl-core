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
[residuals,W0] =  getresiduals(data.X,T,hmm.train.order,hmm.train.orderGL,hmm.train.timelag,hmm.train.lambdaGL);
hmm=initpriors(data.X,T,hmm,residuals);
hmm=initpost(data.X,T,hmm,residuals,Gamma);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmm] = initpriors(X,T,hmm,residuals)
% define priors

if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; 
else orders = []; end 

ndim = size(X,2);

for k=1:hmm.K,
    if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'full')
        defstateprior(k)=struct('sigma',[],'alpha',[],'delta',[],'Omega',[]);
    else
        defstateprior(k)=struct('sigma',[],'alpha',[],'delta',[]);
    end
    defstateprior(k).sigma = struct('Gam_shape',[],'Gam_rate',[]);
    defstateprior(k).sigma.Gam_shape = 0.1*ones(ndim,ndim); %+ 0.05*eye(ndim);
    defstateprior(k).sigma.Gam_rate = 0.1*ones(ndim,ndim);%  + 0.05*eye(ndim);
    defstateprior(k).alpha = struct('Gam_shape',[],'Gam_rate',[]);
    defstateprior(k).alpha.Gam_shape = 0.1;
    defstateprior(k).alpha.Gam_rate = 0.1*ones(1,length(orders));  
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
if ~isfield(hmm,'state') | ~isfield(hmm.state,'prior'),
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmm] = initpost(X,T,hmm,residuals,Gamma)
% Initialising the posteriors

if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end); 
else orders = []; order = 0; end 

Tres = sum(T) - length(T)*order;
N = length(T);
ndim = size(X,2);
K = hmm.K;

% initial random values for the states - multinomial
gammasum = sum(Gamma);

XX = []; Y = [];
for in=1:N
    t0 = sum(T(1:in-1)); s0 = t0 - order*(in-1);
    XX0 = zeros(T(in)-order,length(orders)*ndim);
    for i=1:length(orders)
        o = orders(i);
        XX0(:,(1:ndim) + (i-1)*ndim) = X(t0+order-o+1:t0+T(in)-o,:);
    end;
    XX = [XX; XX0];
    Y = [Y; residuals(s0+1:s0+T(in)-order,:)];
end

% W
for k=1:K,
    hmm.state(k).W = struct('Mu_W',zeros(ndim*length(orders),ndim),'S_W',[]);
    if strcmp(hmm.train.covtype,'uniquediag') || strcmp(hmm.train.covtype,'diag')
        if order==0
            hmm.state(k).W.S_W = zeros(0,0);
        else
            hmm.state(k).W.S_W = zeros(ndim,ndim*length(orders),ndim*length(orders));
            %hmm.state(k).W.S_W(1,:,:) = inv( (XX' .* repmat(Gamma(:,k)',ndim*order,1)) * XX)   ;
            hmm.state(k).W.S_W(1,:,:) = inv( ((XX' .* repmat(Gamma(:,k)',ndim*length(orders),1)) * XX) + 0.01*eye(ndim*length(orders))  )   ;
            for n=2:ndim
                hmm.state(k).W.S_W(n,:,:) = hmm.state(k).W.S_W(1,:,:);
            end;
            hmm.state(k).W.Mu_W = (( permute(hmm.state(k).W.S_W(1,:,:),[2 3 1]) * XX') .* repmat(Gamma(:,k)',ndim*length(orders),1)) * residuals;
        end
    else
        if order==0
            hmm.state(k).W.S_W = zeros(0,0);
        else
            XXGXX = (XX' .* repmat(Gamma(:,k)',ndim*length(orders),1)) * XX;
            hmm.state(k).W.S_W = inv( kron(XXGXX,eye(ndim)) );
            hmm.state(k).W.Mu_W = (( XXGXX \ XX' ) .* repmat(Gamma(:,k)',ndim*length(orders),1)) * residuals;
        end
    end
end;

% Omega
switch hmm.train.covtype
    case 'uniquediag'
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;
        hmm.Omega.Gam_rate = hmm.prior.Omega.Gam_rate;
        for k=1:hmm.K
            e = residuals - XX * hmm.state(k).W.Mu_W;
            hmm.Omega.Gam_rate = hmm.Omega.Gam_rate + 0.5 * sum( repmat(Gamma(:,k),1,ndim) .* e.^2 );
        end
    case 'uniquefull'
        hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres;
        hmm.Omega.Gam_rate = hmm.prior.Omega.Gam_rate;
        for k=1:hmm.K
            e = residuals - XX * hmm.state(k).W.Mu_W;
            hmm.Omega.Gam_rate = hmm.Omega.Gam_rate + (e' .* repmat(Gamma(:,k)',ndim,1)) * e;
        end
        hmm.Omega.Gam_irate = inv(hmm.Omega.Gam_rate);
    case 'diag'
        for k=1:hmm.K
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k) / 2;
            e = residuals - XX * hmm.state(k).W.Mu_W;
            hmm.state(k).Omega.Gam_rate = hmm.state(k).prior.Omega.Gam_rate + 0.5* sum( repmat(Gamma(:,k),1,ndim) .* e.^2 );
        end
    case 'full'
        for k=1:hmm.K
            hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k);
            e = residuals - XX * hmm.state(k).W.Mu_W;
            hmm.state(k).Omega.Gam_rate =  hmm.state(k).prior.Omega.Gam_rate + (e' .* repmat(Gamma(:,k)',ndim,1)) * e;
            hmm.state(k).Omega.Gam_irate = inv(hmm.state(k).Omega.Gam_rate);
        end
end

% sigma
for k=1:K,
    %shape
    hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape;
    hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    hmm.state(k).sigma.Gam_shape = hmm.state(k).prior.sigma.Gam_shape + order;
    for n=1:ndim
        hmm.state(k).sigma.Gam_shape(n,n) = hmm.state(k).prior.sigma.Gam_shape(n,n) + 0.5*order;
    end;
    %rate
    hmm.state(k).sigma.Gam_rate = hmm.state(k).prior.sigma.Gam_rate;
    if order>0
        % mu
        for n1=1:ndim-1
            for n2=n1+1:ndim
                index = n1 + (0:length(orders)-1)*ndim;
                hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                    0.5 * (hmm.state(k).W.Mu_W(index,n2)' * ...
                    ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .*  hmm.state(k).W.Mu_W(index,n2))  );
                index = n2 + (0:length(orders)-1)*ndim;
                hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                    0.5 * (hmm.state(k).W.Mu_W(index,n1)' * ...
                    ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n1)) );
            end
        end
        for n=1:ndim
            index = n + (0:length(orders)-1)*ndim;
            hmm.state(k).sigma.Gam_rate(n,n) = hmm.state(k).sigma.Gam_rate(n,n) + ...
                0.5 * (hmm.state(k).W.Mu_W(index,n)' * ...
                ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n)));
        end;
        % S
        if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
            for n1=1:ndim-1
                for n2=n1+1:ndim
                    index = (0:length(orders)-1) * ndim^2 + (n1-1) * ndim + n2;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * sum( (hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index) ));
                    index = (0:length(orders)-1) * ndim^2 + (n2-1) * ndim + n1;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * sum( (hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index) ));
                    hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n1,n2);
                end
            end
            for n=1:ndim
                index = (0:length(orders)-1) * ndim^2 + (n-1) * ndim + n;
                hmm.state(k).sigma.Gam_rate(n,n) = hmm.state(k).sigma.Gam_rate(n,n) + ...
                    0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index)));
            end
        else
            for n1=1:ndim-1
                for n2=n1+1:ndim
                    index = n1 + (0:length(orders)-1)*ndim;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n2,index,index),[2 3 1])   )) ;
                    index = n2 + (0:length(orders)-1)*ndim;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n1,index,index),[2 3 1])   )) ;
                    hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n1,n2);
                end
            end
            for n=1:ndim
                index = n + (0:length(orders)-1)*ndim;
                hmm.state(k).sigma.Gam_rate(n,n) = hmm.state(k).sigma.Gam_rate(n,n) + ...
                    0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(permute(hmm.state(k).W.S_W(n,index,index),[2 3 1]) ));
            end
        end;
    end;
end;

% alpha
for k=1:K,
    hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape + 0.5*ndim^2;
    hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
        for n=1:ndim,
            for i=1:length(orders)
                index = (i-1)*ndim+n;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( hmm.state(k).W.Mu_W(index,:) .* (hmm.state(k).sigma.Gam_shape(n,:) ./ hmm.state(k).sigma.Gam_rate(n,:) ) ) * ...
                    hmm.state(k).W.Mu_W(index,:)';
                index =  (i-1) * ndim^2 + (n-1) * ndim + (1:ndim);
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * sum(diag(hmm.state(k).W.S_W(index,index)) .* (hmm.state(k).sigma.Gam_shape(:,n) ./ hmm.state(k).sigma.Gam_rate(:,n) ) );
            end;
        end;
    else
        for n=1:ndim,
            for i=1:length(orders),
                index = (i-1)*ndim+n;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( (hmm.state(k).W.Mu_W(index,:) .* (hmm.state(k).sigma.Gam_shape(n,:) ./ hmm.state(k).sigma.Gam_rate(n,:)) ) * ...
                    hmm.state(k).W.Mu_W(index,:)' + sum( (hmm.state(k).sigma.Gam_shape(n,:) ./ hmm.state(k).sigma.Gam_rate(n,:)) .* ...
                    hmm.state(k).W.S_W(:,index,index)'));
            end;
        end;
    end
end;

