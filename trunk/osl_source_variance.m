function V = osl_source_variance(D)
% Computes the temporal variance of data in MEEG object
% V = osl_source_variance(D)
currentMontage = montage(D,'getindex');

V = zeros(D.nchannels,D.ntrials);

D = D.montage('switch');
C = osl_cov(D);

if currentMontage ~= 0
    tra = D.montage('getmontage',currentMontage).tra;
else
    tra = eye(size(C));
end

ft_progress('init', 'etf');
for trl = 1:D.ntrials
    ft_progress(trl/D.ntrials, 'Processing trial %d of %d', trl, D.ntrials);
    V(:,trl) = diag(tra * C(:,:,trl) * tra');
end
ft_progress('close')


end

