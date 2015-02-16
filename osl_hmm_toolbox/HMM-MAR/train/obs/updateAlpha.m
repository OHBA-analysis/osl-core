function hmm = updateAlpha(hmm,orders)

K = length(hmm.state); ndim = size(hmm.state(1).W.Mu_W,2);
S = hmm.train.S;

for k=1:K,
    hmm.state(k).alpha.Gam_shape = hmm.state(k).prior.alpha.Gam_shape;
    hmm.state(k).alpha.Gam_rate = hmm.state(k).prior.alpha.Gam_rate;
    if strcmp(hmm.train.covtype,'full') || strcmp(hmm.train.covtype,'uniquefull')
        for n=1:ndim,
            indr = S(n,:)==1;
            for i=1:length(orders),
                index = (i-1)*ndim+n + ~hmm.train.zeromean;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( hmm.state(k).W.Mu_W(index,indr) .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr) ) ) * ...
                    hmm.state(k).W.Mu_W(index,indr)';
                index =  (i-1) * ndim^2 + (n-1) * ndim + (1:ndim) + (~hmm.train.zeromean)*ndim; index = index(indr);
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * sum(diag(hmm.state(k).W.S_W(index,index))' .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr) ) );
            end;
            hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5 * sum(indr);
        end;
    else
        for n=1:ndim,
            indr = S(n,:)==1;
            for i=1:length(orders),
                index = (i-1)*ndim+n + ~hmm.train.zeromean;
                hmm.state(k).alpha.Gam_rate(i) = hmm.state(k).alpha.Gam_rate(i) + ...
                    0.5 * ( (hmm.state(k).W.Mu_W(index,indr) .* (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr)) ) * ...
                    hmm.state(k).W.Mu_W(index,indr)' + sum( (hmm.state(k).sigma.Gam_shape(n,indr) ./ hmm.state(k).sigma.Gam_rate(n,indr)) .* ...
                    hmm.state(k).W.S_W(indr,index,index)'));
            end;
            hmm.state(k).alpha.Gam_shape = hmm.state(k).alpha.Gam_shape + 0.5 * sum(indr);
        end;
    end
end;

end