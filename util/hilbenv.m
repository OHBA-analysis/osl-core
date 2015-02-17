function [env,t_env] = hilbenv(dat,t,winsize,downsample)
% HILBENV computes the Hilbert envelope of an [N x tpts] matrix, and
% optionally applies a moving average filter.
%
% [env,t_env] = hilbenv(dat,t,winsize,downsample)
% 
% REQUIRED INPUTS:
%
% dat       - An [N x tpts] data matrix.
%
% OPTIONAL INPUTS:
%
% t             - [1 x tpts] vector of time points.
% winsize       - Size of the moving average filter in samples. For no
%                 averaging, specify winsize = 1 or [].
% downsample    - downsamples the output [0/1] (default = 1)
%
% Adam Baker 2014


if ~exist('t','var') || isempty(t)
    t = 1:size(dat,2);
end

if ~exist('winsize','var') || isempty(winsize)
    winsize = 1;
end

if ~exist('downsample','var') || isempty(downsample)
    downsample = 1;
end


ds_fac = ceil(winsize/4);
if ds_fac < 1
    ds_fac = 1;
    t_env = t;
else
    t_env = t(ceil(winsize/2):end-floor(winsize/2));
    t_env = t_env(ds_fac:ds_fac:end);
end


% Compute Hilbert transform:
env = transpose(abs(hilbert(dat')));

% Apply moving average:
if winsize > 1
    
    if 0
        % old code for doing moving averaging and downsampling
        overlap = 0.7500;
        clear env2;
        for vox=1:size(env,1),
            [tmp,t_avg] = osl_movavg(env(vox,:),t,winsize,overlap,0);
            env2(vox,:)=tmp(~isnan(t_avg));
            t_avg=t_avg(~isnan(t_avg));
        end;
        env=env2;
        t_env=t_avg;
    
    else
    
        env = fftfilt(ones(1,winsize),env.', nextpow2(length(env)+1))./winsize;

        % Remove edge effects
        env = env(winsize:end,:);
        env = env(1:end-rem(size(env,1),ds_fac),:);

        t = t(winsize:end);
        t = t(1:end-rem(size(t,1),ds_fac));

        % Downsample the envelope
        if downsample
            env = resample(env,1,ds_fac)';
            env(env<0) = 0;
        else
            t_env = t;
        end
    end;
    
end

end