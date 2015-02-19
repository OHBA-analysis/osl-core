function Denv = osl_hilbenv(S)
% Computes the Hilbert envelope of MEEG data
% Dnew = osl_hilbenv(S)
%
% S.D       - MEEG object
% S.winsize - window size (seconds)
% S.prefix  - filename prefix for new MEEG object
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

S.prefix = ft_getopt(S,'prefix','h');

S.winsize = ft_getopt(S,'winsize');
if isempty(S.winsize)
    S.winsize = 0; % No smoothing
end
  
% Set up new MEEG object to hold downsampled envelope
[~,t_env] = hilbenv(D.time,D.time,round(S.winsize*D.fsample));


Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,S.prefix),[D.nchannels,length(t_env),D.ntrials]);
Denv = timeonset(Denv,t_env(1));

fsampleNew = 1./(diff(t_env([1,end]))/length(t_env));
Denv = fsample(Denv,fsampleNew);

% Loop over blocks and voxels and apply envelope averaging
blks = osl_memblocks(D,1);

ft_progress('init','etf')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    for trl = 1:D.ntrials
        dat_blk = D(blks(iblk,1):blks(iblk,2),:,trl);
        env = hilbenv(dat_blk,D.time,round(S.winsize*D.fsample));   
        Denv(blks(iblk,1):blks(iblk,2),:,trl) = env;
    end
end

ft_progress('close');

Denv.save;

end




