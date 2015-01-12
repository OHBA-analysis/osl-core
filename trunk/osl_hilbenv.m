function Denv = osl_hilbenv(S)
% Computes the Hilbert envelope of MEEG data
% Dnew = osl_hilbenv(S)
%
% S.D       - MEEG object
% S.winsize - window size (samples)
%
% Adam Baker 2014

% Check SPM File Specification:
try
    if isa(S.D,'meeg')
        D = S.D;
    else
        S.D = char(S.D);
        [pathstr,filestr] = fileparts(S.D);
        S.D = fullfile(pathstr,[filestr '.mat']); % force .mat suffix
        D = spm_eeg_load(S.D);
    end;
    D.check;
catch
    error('SPM file specification not recognised or incorrect');
end

S.winsize = ft_getopt(S,'winsize');
if isempty(S.winsize)
    S.winsize = 0; % No smoothing
end

trl = 1; % TODO - fix for trialwise data
    
% Set up new MEEG object to hold downsampled envelope
ds_fac = ceil(S.winsize/2);
if ds_fac < 1
    ds_fac = 1;
    t_env = D.time;
else
    t_env = D.time(ceil(S.winsize/2):end-floor(S.winsize/2));
    t_env = t_env(ds_fac:ds_fac:end);
end


Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,'h'),[D.nchannels,length(t_env),D.ntrials]);
Denv = timeonset(Denv,t_env(1));
Denv = fsample(Denv,Denv.fsample/ds_fac);

% Loop over blocks and voxels and apply envelope averaging
blks = osl_memblocks(D,1);

ft_progress('init','eta')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    % Compute Hilbert transform:
    dat_blk = transpose(abs(hilbert(D(blks(iblk,1):blks(iblk,2),:,trl)')));
  
    % Apply moving average:
    if S.winsize > 0
        dat_blk = fftfilt(repmat(ones(S.winsize,1),1,size(dat_blk,1)),dat_blk');
        dat_blk = dat_blk./S.winsize;
        dat_blk = dat_blk(S.winsize:end,:);
        dat_blk = dat_blk(1:end-rem(size(dat_blk,1),ds_fac),:);
        
        % Downsample
        dat_blk = resample(dat_blk,1,ds_fac)';
        dat_blk(dat_blk<0) = 0;
    end
    
    Denv(blks(iblk,1):blks(iblk,2),:) = dat_blk;
end
ft_progress('close');

Denv.save;


end




