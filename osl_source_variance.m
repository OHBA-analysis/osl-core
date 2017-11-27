function [V, mean_V] = osl_source_variance(D)
% Computes the temporal variance of data in MEEG object
% V = osl_source_variance(D)
% [V, mean_V] = osl_source_variance(D)

currentMontage = montage(D,'getindex');

V = zeros(D.nchannels,D.ntrials);

if nargout>1
   mean_V=zeros(D.nchannels,1); 
   num_good_trials=0;
end

D = D.montage('switch');
C = osl_cov(D);

if currentMontage ~= 0
    tra = D.montage('getmontage',currentMontage).tra;
else
    tra = eye(size(C(:,:,1)));
end

if D.ntrials > 1, ft_progress('init', 'etf'); end
for trl = 1:D.ntrials
    if D.ntrials > 1
        ft_progress(trl/D.ntrials, 'Processing trial %d of %d', trl, D.ntrials); 
    end
    V(:,trl) = diag(tra * C(:,:,trl) * tra');
    
    if nargout>1
        % compute source variance averaged over good trials
        if sum(V(:,trl))>0 % check if trial was good
            mean_V=mean_V+V(:,trl);
            num_good_trials=num_good_trials+1;
        end
    end 
end

if D.ntrials > 1, ft_progress('close'); end

if nargout>1
    % compute source variance averaged over trials        
    mean_V = mean_V/num_good_trials;
end

end

