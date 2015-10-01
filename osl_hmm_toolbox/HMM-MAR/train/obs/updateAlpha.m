function hmm = updateAlpha(hmm)

K = length(hmm.state); ndim = hmm.train.ndim;

for k=1:K,
    setstateoptions;
    if isempty(orders) || ~isempty(train.prior), continue; end
    hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape;
    hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    if train.uniqueAR 
        hmm.state(k).alpha.Gam_rate = hmm.state(k).alpha.Gam_rate + ...
            0.5 * ( hmm.state(k).W.Mu_W(1+(~train.zeromean) : end).^2 )' + ...
            0.5 * diag(hmm.state(k).W.S_W)' ;
        hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5;
        
    elseif ndim==1
        hmm.state(k).alpha.Gam_rate = hmm.state(k).alpha.Gam_rate + ...
            0.5 * ( hmm.state(k).W.Mu_W(1+(~train.zeromean) : end).^2 )' + ...
            0.5 * diag(permute(hmm.state(k).W.S_W,[2 3 1]))' ;
        hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5;
   
    elseif strcmp(train.covtype,'full') || strcmp(train.covtype,'uniquefull')
        for n=1:ndim,
            indr = S(n,:)==1;
            for i=1:length(orders),
                index = (i-1)*ndim+n + ~train.zeromean;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( hmm.state(k).W.Mu_W(index,indr) .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr) ) ) * ...
                    hmm.state(k).W.Mu_W(index,indr)';
                index =  (i-1) * ndim^2 + (n-1) * ndim + (1:ndim) + (~train.zeromean)*ndim; index = index(indr);
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * sum(diag(hmm.state(k).W.S_W(index,index))' .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr) ) );
            end
            hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5 * sum(indr);
        end
        
    else
        for n=1:ndim,
            indr = S(n,:)==1;
            for i=1:length(orders),
                index = (i-1)*ndim+n + ~train.zeromean;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( (hmm.state(k).W.Mu_W(index,indr) .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr)) ) * ...
                    hmm.state(k).W.Mu_W(index,indr)' + sum( (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr)) .* ...
                    hmm.state(k).W.S_W(indr,index,index)'));
            end;
            hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5 * sum(indr);
        end
    end
end

end