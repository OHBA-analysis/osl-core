function [env,t_env] = hilbenv(dat,t,winsize,badsamples)
% HILBENV computes the Hilbert envelope of an [N x tpts] matrix, and
% optionally applies a moving average filter.
%
% [env,t_env] = hilbenv(dat,t,winsize)
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
% badsamples    - list of bad samples or []
%
% Adam Baker 2014


if ~exist('t','var') || isempty(t)
    t = 1:size(dat,2);
end

if ~exist('winsize','var') || isempty(winsize)
    winsize = 1;
end

ds_fac = ceil(winsize/4);
if ds_fac < 1
    t_env = t;
else
    t_env = t(ceil(winsize/2):end-floor(winsize/2));
    t_env = t_env(ds_fac:ds_fac:end);
end


% Compute Hilbert transform:
env = transpose(abs(hilbert(dat')));

% Replace badsamples with NaN:
env(:,badsamples) = nan;

% Apply moving average:
if winsize > 1
    
    overlap = 0.75;
    
    [~,t_avg] = movavg(env(1,:),t,winsize,overlap,0);
    t_avg = t_avg(~isnan(t_avg));
    
    env2 = zeros(size(env,1),length(t_avg));
    
    for vox = 1:size(env,1),
        [tmp,~] = movavg(env(vox,:),t,winsize,overlap,0);
        env2(vox,:) = tmp(~isnan(tmp));
    end
    
    env   = env2;
    t_env = t_avg;
    
end

end