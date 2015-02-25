function [data_avg,t_avg] = osl_movavg(data,t,winsize,overlap,resamp,robust)
% Moving window averaging
% [data_avg,t_avg] = OSL_MOVAVG(data,t,winsize,overlap)
% -----------------------------------------------------------------
% data    - data to average
% t       - vector of sample times
% winsize - window to average over (samples)
% overlap - amount of segment overlapping (e.g. 0.75 = 75% overlap)
% resamp  - resample to original data length
% robust  - use robust averaging (% of data to exlude)
% -----------------------------------------------------------------
% AB 2011

use_hanning = 0;
    
if nargin < 4
  overlap = 0.75;
end
if nargin < 5
  resamp = 0;
end
if nargin < 6
  robust = 0;
else
  use_hanning = 0;
end

if isempty(t) 
  t = 1:length(data);
end

if use_hanning
  winsize = winsize *2;
end

data = data(:)'; t = t(:)';

data_pad = [nan(1,winsize),data,nan(1,winsize)];
t_pad    = [nan(1,winsize),t,nan(1,winsize)];

bf       = buffer(data_pad,winsize,round(winsize*overlap),'nodelay');
tbf      = buffer(t_pad,winsize,round(winsize*overlap),'nodelay');

% reject data sections where more than 25% of data is NaN
nans = sum(isnan(bf)) > 0.25*winsize;
bf( :,nans) = nan; % - do this by replacing all with nans
tbf(:,nans) = nan; % - do this by replacing all with nans

if(use_hanning)
hanning_window = repmat(hanning(winsize),1,size(bf,2));
bf = bf.*hanning_window;
end

if robust ~= 0
  data_avg = trimmean(bf,robust);
else
  data_avg = nanmean(bf);
end

t_avg    = nanmean(tbf);

if resamp
  t_start = t_avg(find(~isnan(t_avg),1,'first'));
  t_end   = t_avg(find(~isnan(t_avg),1,'last'));

  data_i = nan(size(t));
  data_i(t>=t_start & t<t_end) = interp1(t_avg(t_avg>=t_start & t_avg<t_end),data_avg(t_avg>=t_start & t_avg<t_end),t(t>=t_start & t<t_end),'linear',nan);
  
  data_avg = data_i;
  data_avg(isnan(data)) = nan;
  t_avg = t;
end

end

