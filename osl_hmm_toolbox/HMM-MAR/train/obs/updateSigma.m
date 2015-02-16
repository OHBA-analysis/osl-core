function hmm = updateSigma(hmm,orders)

K = length(hmm.state); ndim = size(hmm.state(1).W.Mu_W,2);
S = hmm.train.S;

for k=1:K,
    %shape
    if hmm.train.symmetricprior,
        hmm.state(k).sigma.Gam_shape = hmm.state(k).prior.sigma.Gam_shape + length(orders);
        for n=1:ndim
            hmm.state(k).sigma.Gam_shape(n,n) = hmm.state(k).prior.sigma.Gam_shape(n,n) + 0.5*length(orders);
        end
    else
        hmm.state(k).sigma.Gam_shape = hmm.state(k).prior.sigma.Gam_shape + 0.5*length(orders);
    end
    %rate
    hmm.state(k).sigma.Gam_rate = hmm.state(k).prior.sigma.Gam_rate;
    % mean(W)
    for n1=1:ndim
        if any(S(n1,:)==1)
            for n2=find(S(n1,:)==1)
                if hmm.train.symmetricprior && n1>n2, 
                    continue; 
                end
                if n1==n2,
                    index = n1 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                    hmm.state(k).sigma.Gam_rate(n1,n1) = hmm.state(k).sigma.Gam_rate(n1,n1) + ...
                        0.5 * (hmm.state(k).W.Mu_W(index,n1)' * ...
                        ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n1)));
                else
                    index = n1 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                    hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                        0.5 * (hmm.state(k).W.Mu_W(index,n2)' * ...
                        ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n2)) );
                    index = n2 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                    h = 0.5 * (hmm.state(k).W.Mu_W(index,n1)' * ...
                        ((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* hmm.state(k).W.Mu_W(index,n1)));
                    if hmm.train.symmetricprior,
                        hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + h;
                    else
                        hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n2,n1) + h;
                    end
                end
            end
        end
    end
    % cov(W)
    if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
        for n1=1:ndim
            if any(S(n1,:)==1)
                for n2=find(S(n1,:)==1)
                    if hmm.train.symmetricprior && n1>n2,
                        continue;
                    end
                    if n1==n2,
                        index = (0:length(orders)-1) * ndim^2 + (n1-1) * ndim + n2 + (~hmm.train.zeromean)*ndim;
                        hmm.state(k).sigma.Gam_rate(n1,n1) = hmm.state(k).sigma.Gam_rate(n1,n1) + ...
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index)));
                    else
                        index = (0:length(orders)-1) * ndim^2 + (n1-1) * ndim + n2 + (~hmm.train.zeromean)*ndim;
                        hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index) ));
                        index = (0:length(orders)-1) * ndim^2 + (n2-1) * ndim + n1 + (~hmm.train.zeromean)*ndim;
                        h = 0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag(hmm.state(k).W.S_W(index,index) ));
                        if hmm.train.symmetricprior,
                            hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + h;
                        else
                            hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n2,n1) + h;
                        end
                    end
                end
            end
        end

    else
        for n1=1:ndim
            if any(S(n1,:)==1)
                for n2=find(S(n1,:)==1)
                    if hmm.train.symmetricprior && n1>n2,
                        continue;
                    end
                    if n1==n2,
                        index = n1 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                        hmm.state(k).sigma.Gam_rate(n1,n1) = hmm.state(k).sigma.Gam_rate(n1,n1) + ...
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* diag( permute(hmm.state(k).W.S_W(n1,index,index),[2 3 1]) ));
                    else
                        index = n1 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                        hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + ...
                            0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* ...
                            diag( permute(hmm.state(k).W.S_W(n2,index,index),[2 3 1]) )) ;
                        index = n2 + (0:length(orders)-1)*ndim + ~hmm.train.zeromean;
                        h = 0.5 * sum((hmm.state(k).alpha.Gam_shape ./ hmm.state(k).alpha.Gam_rate') .* ...
                            diag( permute(hmm.state(k).W.S_W(n1,index,index),[2 3 1]) )) ;
                        if hmm.train.symmetricprior,
                            hmm.state(k).sigma.Gam_rate(n1,n2) = hmm.state(k).sigma.Gam_rate(n1,n2) + h;
                        else
                            hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n2,n1) + h;
                        end
                    end
                end
            end
        end
    end;
    if hmm.train.symmetricprior,
        for n1=1:ndim-1
            if any(S(n1,:)==1)
                for n2=find(S(n1,:)==1)
                    if n1==n2, continue; end
                    hmm.state(k).sigma.Gam_rate(n2,n1) = hmm.state(k).sigma.Gam_rate(n1,n2);
                end
            end
        end
    end
end;
end