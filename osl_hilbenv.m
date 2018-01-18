function Denv = osl_hilbenv(S)
% Computes the Hilbert envelope of MEEG data
% Dnew = osl_hilbenv(S)
%
% S.D         - MEEG object
% S.winsize   - window size (seconds)
% S.prefix    - filename prefix for new MEEG object
% S.freqbands - cell array of frequency bands to use [Hz]
%                 i.e. {[f1_low f1_high],[f2_low f2_high]}
% S.logtrans  - apply log transform [0/1] (default 0)
% S.demean    - remove mean from envelope (default 0)
% S.downsample- downsample envelope (default 1)
% S.remove_edge_effects- remove_edge_effects (default 1)

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

S.logtrans = ft_getopt(S,'logtrans',0);
S.demean   = ft_getopt(S,'demean',0);
S.downsample = ft_getopt(S,'downsample',1);
S.remove_edge_effects = ft_getopt(S,'remove_edge_effects',1);

S.winsize = ft_getopt(S,'winsize');
if isempty(S.winsize)
    S.winsize = 0; % No smoothing
end


S.freqbands = ft_getopt(S,'freqbands',{});
if ~iscell(S.freqbands)
    error('S.freqbands must be a cell array');
end
  
% Set up new MEEG object to hold downsampled envelope
[~,t_env] = hilbenv(D.time,D.time,round(S.winsize*D.fsample),S.downsample,S.remove_edge_effects);

if numel(S.freqbands) < 2
    % Create Nchannels x Nsamples x Ntrials object
    Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,S.prefix),[D.nchannels,length(t_env),D.ntrials]);
    Denv = timeonset(Denv,t_env(1));
else
    % Create Nchannels x Nfrequencies x Nsamples x Ntrials TF object
    Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,S.prefix),[D.nchannels,numel(S.freqbands),length(t_env),D.ntrials]);
    Denv = timeonset(Denv,t_env(1));
end

if S.downsample
    fsampleNew = 1./(diff(t_env([1,end]))/length(t_env));
else
    fsampleNew = D.fsample;
end

Denv = fsample(Denv,fsampleNew);

% Loop over blocks and voxels and apply envelope averaging
blks = osl_memblocks(size(D),1,50*2^20);

ft_progress('init','etf')
for iblk = 1:size(blks,1)
    ft_progress(iblk/size(blks,1));

    for trl = 1:D.ntrials
       
        if ~isempty(S.freqbands)
            % Filter and compute envelope for each band
            for f = 1:numel(S.freqbands)
                dat_blk = bandpass(D(blks(iblk,1):blks(iblk,2),:,trl),S.freqbands{f},D.fsample);
                env = hilbenv(dat_blk,D.time,round(S.winsize*D.fsample),S.downsample,S.remove_edge_effects);
                env = permute(env,[1 3 2]);
                if S.logtrans
                    env = log10(env);
                end
                if S.demean
                    env = bsxfun(@minus,env,mean(env,2));
                end
                Denv(blks(iblk,1):blks(iblk,2),f,:,trl) = env;
            end
        else
            % Compute wideband envelope
            dat_blk = D(blks(iblk,1):blks(iblk,2),:,trl); 
            env = hilbenv(dat_blk,D.time,round(S.winsize*D.fsample),S.downsample,S.remove_edge_effects);
            if S.logtrans
                env = log10(env);
            end
            if S.demean
                env = bsxfun(@minus,env,mean(env,2));
            end
            Denv(blks(iblk,1):blks(iblk,2),:,trl) = env;
        end
        
    end
end

ft_progress('close');

Denv.save;

end




