function [hmm] = obsupdate (X,T,Gamma,hmm,residuals)
%
% Update observation model
%
% INPUT
% X             observations
% T             length of series
% Gamma         p(state given X)
% hmm           hmm data structure
% residuals     in case we train on residuals, the value of those.
%
% OUTPUT
% hmm           estimated HMMMAR model 
%
% Author: Diego Vidaurre, OHBA, University of Oxford

if hmm.train.order > 0, orders = 1:hmm.train.timelag:hmm.train.order; order = orders(end); 
else orders = []; order = 0; end 

Tres = sum(T) - length(T)*order;
N = length(T);
ndim = size(X,2);
K=hmm.K;

obs_tol = 0.00005;
obs_maxit = 2; %20;
mean_change = Inf;
obs_it = 1;

% Some stuff that will be later used
M_Gamma2 = gamma_square(Gamma);
Gamma2 = zeros(K,1);
for k=1:K,
    Gamma2(k) = sum(M_Gamma2(:,k));
end;
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

XXGXX = zeros(length(orders)*ndim,length(orders)*ndim,K);
for k=1:K
    XXGXX(:,:,k) = (XX' .* repmat(Gamma(:,k)',ndim*length(orders),1)) * XX;
end

if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
    mlW = zeros(length(orders)*ndim*ndim,K);
    if order>0
        for k=1:K
            tmp = (( XXGXX(:,:,k) \ XX') .* repmat(Gamma(:,k)',ndim*length(orders),1) * residuals)';
            mlW(:,k) = tmp(:);
        end
    end
end


while mean_change>obs_tol && obs_it<=obs_maxit,
    
    last_state = hmm.state;
    
    %%% W
    if 1
        if order>0
            if strcmp(hmm.train.covtype,'diag') || strcmp(hmm.train.covtype,'uniquediag')
                for k=1:K,
                    alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate), ndim, 1);
                    for n=1:ndim
                        if strcmp(hmm.train.covtype,'diag')
                            omega = hmm.state(k).Omega;
                        else
                            omega = hmm.Omega;
                        end
                        regterm = diag( repmat(hmm.state(k).sigma.Gam_shape(n,:) ./  hmm.state(k).sigma.Gam_rate(n,:), 1, length(orders)) .* alphaterm(:)' );
                        iS =  regterm + (omega.Gam_shape / omega.Gam_rate(n)) * XXGXX(:,:,k);
                        S = inv(iS);
                        hmm.state(k).W.S_W(n,:,:) = S;
                        hmm.state(k).W.Mu_W(:,n) = ((S * (omega.Gam_shape / omega.Gam_rate(n)) * XX') .* repmat(Gamma(:,k)',ndim*length(orders),1)) * residuals(:,n);
                    end;
                end;
            else
                for k=1:K
                    alphaterm = repmat( (hmm.state(k).alpha.Gam_shape ./  hmm.state(k).alpha.Gam_rate), ndim^2, 1);
                    alphaterm = alphaterm(:);
                    sigmaterm = (hmm.state(k).sigma.Gam_shape ./ hmm.state(k).sigma.Gam_rate);
                    sigmaterm = repmat(sigmaterm(:), length(orders), 1);
                    regterm = diag(alphaterm .* sigmaterm);
                    if strcmp(hmm.train.covtype,'full')
                        omega = hmm.state(k).Omega;
                    else
                        omega = hmm.Omega;
                    end
                    prec = omega.Gam_shape * omega.Gam_irate;
                    gram = kron(XXGXX(:,:,k), prec);
                    hmm.state(k).W.S_W = inv(regterm + gram);
                    muW = hmm.state(k).W.S_W * gram * mlW(:,k);
                    hmm.state(k).W.Mu_W = reshape(muW,ndim,ndim*length(orders))';
                end
            end;
        end
    end
    
    XW = zeros(K,Tres,ndim);
    for k=1:K
        XW(k,:,:) = XX * hmm.state(k).W.Mu_W;
    end;
    
    %%% Omega
    if 1
        switch hmm.train.covtype,
            case 'uniquediag'
                hmm.Omega.Gam_rate = hmm.prior.Omega.Gam_rate;
                for k=1:K
                    XW2d = permute(XW(k,:,:),[2 3 1]);
                    e = (residuals - XW2d).^2;
                    swx2 = zeros(Tres,ndim);
                    if order>0
                        for n=1:ndim
                            tmp = XX * permute(hmm.state(k).W.S_W(n,:,:),[2 3 1]);
                            swx2(:,n) = sum(tmp .* XX,2);
                        end;
                    end
                    hmm.Omega.Gam_rate = hmm.Omega.Gam_rate +  0.5* sum( repmat(Gamma(:,k),1,ndim) .* (e + swx2) );
                end;
                hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres / 2;
                
            case 'diag'
                for k=1:K
                    XW2d = permute(XW(k,:,:),[2 3 1]);
                    e = (residuals - XW2d).^2;
                    swx2 = zeros(Tres,ndim);
                    if order>0
                        for n=1:ndim
                            tmp = XX * permute(hmm.state(k).W.S_W(n,:,:),[2 3 1]);
                            swx2(:,n) = sum(tmp .* XX,2);
                        end;
                    end
                    hmm.state(k).Omega.Gam_rate = hmm.state(k).prior.Omega.Gam_rate + ...
                        0.5* sum( repmat(Gamma(:,k),1,ndim) .* (e + swx2) );
                    hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k) / 2;
                end
                
            case 'uniquefull'
                hmm.Omega.Gam_rate = hmm.prior.Omega.Gam_rate;
                for k=1:K
                    XW2d = permute(XW(k,:,:),[2 3 1]);
                    e = (residuals - XW2d);
                    e = (e' .* repmat(Gamma(:,k)',ndim,1)) * e;
                    swx2 =  zeros(ndim,ndim);
                    if order>0
                        for n1=1:ndim
                            for n2=n1:ndim
                                index1 = (0:length(orders)*ndim-1) * ndim + n1;
                                index2 = (0:length(orders)*ndim-1) * ndim + n2;
                                swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX(:,:,k)));
                                swx2(n2,n1) = swx2(n1,n2);
                            end
                        end
                    end
                    hmm.Omega.Gam_rate = hmm.Omega.Gam_rate +  (e + swx2);
                end
                hmm.Omega.Gam_irate = inv(hmm.Omega.Gam_rate);
                hmm.Omega.Gam_shape = hmm.prior.Omega.Gam_shape + Tres;
                
            case 'full'
                for k=1:K
                    XW2d = permute(XW(k,:,:),[2 3 1]);
                    e = (residuals - XW2d);
                    e = (e' .* repmat(Gamma(:,k)',ndim,1)) * e;
                    swx2 =  zeros(ndim,ndim);
                    if order>0
                        for n1=1:ndim
                            for n2=n1:ndim
                                index1 = (0:length(orders)*ndim-1) * ndim + n1;
                                index2 = (0:length(orders)*ndim-1) * ndim + n2;
                                swx2(n1,n2) = sum(sum(hmm.state(k).W.S_W(index1,index2) .* XXGXX(:,:,k)));
                                swx2(n2,n1) = swx2(n1,n2);
                            end
                        end
                    end
                    hmm.state(k).Omega.Gam_rate = hmm.state(k).prior.Omega.Gam_rate + (e + swx2);
                    hmm.state(k).Omega.Gam_shape = hmm.state(k).prior.Omega.Gam_shape + gammasum(k);
                    hmm.state(k).Omega.Gam_irate = inv(hmm.state(k).Omega.Gam_rate);
                end
        end;
    end;
    
    %%% sigma - channel x channel coefficients (symmetric)
    if order>=1
        for k=1:K,
            %shape
            hmm.state(k).sigma.Gam_shape = hmm.state(k).prior.sigma.Gam_shape + length(orders);
            for n=1:ndim
                hmm.state(k).sigma.Gam_shape(n,n) = hmm.state(k).prior.sigma.Gam_shape(n,n) + 0.5*length(orders);
            end;
            %rate
            hmm.state(k).sigma.Gam_rate = hmm.state(k).prior.sigma.Gam_rate;
            % mu
            for n1=1:ndim-1
                for n2=n1+1:ndim
                    index = n1 + (0:length(orders)-1)*ndim;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * (hmm.state(k).W.Mu_W(index,n2)' * ...
                        ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n2)) );
                    index = n2 + (0:length(orders)-1)*ndim;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * (hmm.state(k).W.Mu_W(index,n1)' * ...
                        ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n1)));
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
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index) ));
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
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n2,index,index),[2 3 1]) )) ;
                        index = n2 + (0:length(orders)-1)*ndim;
                        hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n1,index,index),[2 3 1]) )) ;
                        hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n1,n2);
                    end
                end
                for n=1:ndim
                    index = n + (0:length(orders)-1)*ndim;
                    hmm.state(k).sigma.Gam_rate(n,n) = hmm.state(k).sigma.Gam_rate(n,n) + ...
                        0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n,index,index),[2 3 1]) ));
                end
            end;
        end;
    end;
    
    
    %%% alpha: one per order
    if order>=1
        for k=1:K,
            hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape + 0.5*ndim^2;
            hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
            if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
                for n=1:ndim,
                    for i=1:length(orders),
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
    end
    
    %%% termination conditions
    obs_it = obs_it + 1;
    mean_changew = 0;
    for k=1:K
        mean_changew = mean_changew + sum(sum(abs(last_state(k).W.Mu_W - hmm.state(k).W.Mu_W))) / length(orders) / ndim^2 / K;
    end;
    mean_change = mean_changew;
end;

