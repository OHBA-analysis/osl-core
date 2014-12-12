function Denv = osl_hilbenv(S)
% Computes the Hilbert envelope of MEEG data
% Dnew = osl_hilbenv(S)
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
    S.winsize = D.fsample; % 1 second
end

S.overlap = ft_getopt(S,'overlap');
if isempty(S.overlap)
    S.overlap = 0.75;
end

trl = 1; % TODO - fix for trialwise data

Mem_max = 200*2^20;
Mem_chan = 8*numel(D(1,:,trl));
blk_size = floor(Mem_max./Mem_chan);
    
% Set up new MEEG object to hold downsampled envelope
ds_fac = ceil(S.winsize/2);

t_env = D.time(floor(S.winsize/2):end-ceil(S.winsize/2));
t_env = t_env(ds_fac:ds_fac:end);

Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,'h'),[D.nchannels,length(t_env),D.ntrials]);
Denv = timeonset(Denv,t_env(1));
Denv = fsample(Denv,Denv.fsample/ds_fac);

% Loop over blocks and voxels and apply envelope averaging
blks = 1:blk_size:D.nchannels;
blks = unique([blks D.nchannels]);
blks = [blks(1:end-1); blks(2:end)]';

ft_progress('init','eta')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    % Load in data for this block and rearrange into complex form
    dat_blk = transpose(abs(hilbert(D(blks(iblk,1):blks(iblk,2),:,trl)')));
    
    dat_blk = fftfilt(repmat(ones(S.winsize,1),1,size(dat_blk,1)),dat_blk');
    dat_blk = dat_blk./S.winsize;
    dat_blk = dat_blk(S.winsize:end,:);
    dat_blk = dat_blk(1:end-rem(size(dat_blk,1),ds_fac),:);
    
    % Downsample
    dat_blk = resample(dat_blk,1,ds_fac)';
    
    Denv(blks(iblk,1):blks(iblk,2),:) = dat_blk;
end
ft_progress('close');

Denv.save;


end




