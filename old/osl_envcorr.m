function C = osl_envcorr(S)
% Computes the envelope correlation between a seed coordinate and multiple
% target coordinates.
%
% C = osl_envcorr(S)
% 
% D             - MEEG object containing the data
% seedind       - index into channel/voxel in D for the seed
% targetind     - indices into channel(s)/voxel(s) for the target
% winsize       - moving average window size to use (samples)
% orthogonalise - apply orthogonalisation to remove signal leakage [0/1]
%
% Adam Baker 2014

D = spm_eeg_load(S.D);

trl = 1; % TODO - fix for trialwise data

samples2use = find(goodsamples(D));

% Check 
try prcorr(1,1);
    useprcorr = 1;
catch
    useprcorr = 0;
end


% Compute Hilbert transform for seed:
seed_data = D(S.seedind,samples2use,trl);
seed_env  = hilbenv(seed_data,D.time,S.winsize,1);


% Read data from file in blocks
blks = osl_memblocks(size(D),1);

C = zeros(length(S.targetind),D.ntrials);

ft_progress('init','eta')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    dat_blk = D(S.targetind(blks(iblk,1):blks(iblk,2)),samples2use,trl);

    
    % Apply orthogonalisation
    if S.orthogonalise
        for i = 1:size(dat_blk,1)
            dat_blk(i,:) = dat_blk(i,:)' - seed_data(:)*(seed_data(:)\dat_blk(i,:)');
        end
    end
    
    % Compute Hilbert transform:
    dat_blk = hilbenv(dat_blk,D.time,S.winsize,1);
    
    % Correlate
    for i = 1:size(dat_blk,1)
        if useprcorr
            C(blks(iblk,1)+i-1,trl) = prcorr(seed_env',dat_blk(i,:)');
        else
            C(blks(iblk,1)+i-1,trl) = corr(seed_env',dat_blk(i,:)');
        end
    end

end
ft_progress('close')

