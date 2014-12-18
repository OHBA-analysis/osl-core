function C = osl_envcorr(D,seedind,targetind,winsize,samples2use)
% Computes the Hilbert envelope of MEEG data
% Dnew = osl_hilbenv(S)
%
% S.D       - MEEG object
% S.winsize - window size (seconds)
%
% Adam Baker 2014

S.orthogonalise = 0;

S.winsize = D.fsample;
ds_fac = ceil(S.winsize/2);

trl = 1; % TODO - fix for trialwise data

if nargin < 5
    samples2use = true(1,D.nsamples);
end

samples2use = samples2use(:)';

samples2use = samples2use & ~osl_bad_sections(D,'logical');
samples2use = find(samples2use);
    
% Check 
try prcorr(1,1);
    useprcorr = 1;
catch
    useprcorr = 0;
end


% Reference channel:

% Compute Hilbert transform:
seed_data = D(seedind,samples2use,trl)';
seed_env = transpose(abs(hilbert(seed_data)));

% Apply moving average:
seed_env = fftfilt(repmat(ones(S.winsize,1),1,size(seed_env,1)),seed_env');
seed_env = seed_env./S.winsize;
seed_env = seed_env(S.winsize:end,:);
seed_env = seed_env(1:end-rem(size(seed_env,1),ds_fac),:);

% Downsample
seed_env = resample(seed_env,1,ds_fac);
seed_env(seed_env<0) = 0;


% Read data from file in blocks
Mem_max = 200*2^20;
Mem_chan = 8*numel(D(1,samples2use,trl));
blk_size = floor(Mem_max./Mem_chan);
blks = 1:blk_size:length(targetind);
blks = unique([blks length(targetind)]);
blks = [blks(1:end-1); blks(2:end)]';

C = zeros(length(targetind),D.ntrials);

ft_progress('init','eta')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    dat_blk = D(targetind(blks(iblk,1):blks(iblk,2)),samples2use,trl)';

    
    % Apply orthogonalisation
    if S.orthogonalise
        for i = 1:size(dat_blk,2)
            dat_blk(:,i) = dat_blk(:,i) - seed_data*(seed_data\dat_blk(:,i));
        end
    end
    % Compute Hilbert transform:
    dat_blk = transpose(abs(hilbert(dat_blk)));
  
    % Apply moving average:
    dat_blk = fftfilt(repmat(ones(S.winsize,1),1,size(dat_blk,1)),dat_blk');
    dat_blk = dat_blk./S.winsize;
    dat_blk = dat_blk(S.winsize:end,:);
    dat_blk = dat_blk(1:end-rem(size(dat_blk,1),ds_fac),:);
    
    % Downsample
    dat_blk = resample(dat_blk,1,ds_fac);
    dat_blk(dat_blk<0) = 0;

    % Correlate
    for i = 1:size(dat_blk,2)
        if useprcorr
            C(blks(iblk,1)+i-1,trl) = prcorr(seed_env,dat_blk(:,i));
        else
            C(blks(iblk,1)+i-1,trl) = corr(seed_env,dat_blk(:,i));
        end
    end

end