function V = glean_variance(D)
% Computes the temporal variance of data in MEEG object
% V = glean_variance(D)
currentMontage = montage(D,'getindex');

V = zeros(D.nchannels,D.ntrials);

D = D.montage('switch');
C = glean_cov(D);

if currentMontage ~= 0
    tra = D.montage('getmontage',currentMontage).tra;
else
    tra = eye(size(C));
end

if D.ntrials > 1, ft_progress('init', 'etf'); end
for trl = 1:D.ntrials
    if D.ntrials > 1
        ft_progress(trl/D.ntrials, 'Processing trial %d of %d', trl, D.ntrials); 
    end
    V(:,trl) = diag(tra * C(:,:,trl) * tra');
end
if D.ntrials > 1, ft_progress('close'); end


end

