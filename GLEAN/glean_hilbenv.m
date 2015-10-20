function Denv = glean_hilbenv(S)
% Optimised Hilbert envelope computation for of MEEG data. 
% Speed improvements come from using resampling (with anti-aliasing) and 
% fast data writing via dataio.m for writing intermediate data to disk
%
% Dnew = GLEAN_HILBENV(S)
%
% REQUIRED INPUTS:
%   S.D           - MEEG object
%
% OPTIONAL INPUTS:
%   S.fsample_new - New sampling rate of the envelope (default original)
%   S.prefix      - filename prefix for new MEEG object (default 'h')
%   S.freqbands   - cell array of frequency bands to use [Hz]
%                     i.e. {[f1_low f1_high],[f2_low f2_high]} 
%                    (default [0 Inf])
%   S.logtrans    - apply log transform [0/1] (default 0)
%   S.demean      - remove mean from envelope (default 0)
%
% OUTPUTS:
%   Denv          - Newly created SPM12 MEEG object containing envelopes
%
% Adam Baker 2015


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

fsample_new = ft_getopt(S,'fsample_new');
if isempty(S.fsample_new)
    fsample_new = D.fsample; % No smoothing
end

if rem(fsample_new,1) ~= 0 
    error('fsample_new must be an integer');
end

if rem(D.fsample,1) ~= 0
    error('D.fsample must be an integer');
end

S.freqbands = ft_getopt(S,'freqbands',{[0 Inf]});
if ~iscell(S.freqbands)
    error('S.freqbands must be a cell array');
end
  
bad_samples = find(all(badsamples(D,':',':',':')));
bad_samples = unique([1 bad_samples D.nsamples]); % ensure edge effects are removed

trl     = 1; % Sort this out!
Nvoxels = D.nchannels;

% Get size of downsampled envelope
t_env = D.time;
t_env(bad_samples) = nan;
t_env = resample(t_env,fsample_new,D.fsample);
t_env = t_env(~isnan(t_env));

% Create Nchannels x Nfrequencies x Nsamples x Ntrials TF object
Denv = clone(montage(D,'switch',0),prefix(D.fnamedat,S.prefix),[Nvoxels,numel(S.freqbands),length(t_env),D.ntrials]);
Denv = timeonset(Denv,t_env(1));
Denv = events(Denv,1,[]); %won't be needing these
Denv = fsample(Denv,fsample_new);

% Make a temporary filename for each frequency band to hold filtered data
for f = 1:numel(S.freqbands)
    
    blks = memblocks(size(D),1);
    
    [~,tempfile] = fileparts(tempname);
    tempfile = fullfile(D.path,[tempfile '.bin']);
    
    disp(['Computing envelopes for band ' num2str(f)])
    ft_progress('init','etf')
    for iblk = 1:size(blks,1)
        
        ft_progress(iblk/size(blks,1));
        
        % Filter:
        if all(isfinite(S.freqbands{f}))
            dat_blk = ft_preproc_bandpassfilter(D(blks(iblk,1):blks(iblk,2),:,trl),D.fsample,S.freqbands{f},5,'but','twopass','reduce');
        else
            dat_blk = D(blks(iblk,1):blks(iblk,2),:,trl);
        end
        
        % Hilbert envelope
        dat_blk = transpose(dat_blk);
        dat_blk = abs(hilbert(dat_blk));
        dat_blk(bad_samples,:) = nan;

        % Downsample envelope
        env = zeros(length(t_env),size(dat_blk,2));
        for vox = 1:size(dat_blk,2)
            tmp = resample(dat_blk(:,vox),fsample_new,D.fsample);
            env(:,vox) = tmp(~isnan(tmp));
            if S.logtrans
                env(:,vox) = log10(env(:,vox));
            end
            if S.demean
                env(:,vox) = env(:,vox) - mean(env(:,vox));
            end
        end
        
        dataio(tempfile,env);
    end
    
    ft_progress('close');
    
    Denv(:,f,:,1) = permute(dataio(tempfile),[2,3,1,4]);
    system(['rm ' tempfile]);
    
end

Denv.save;

end






