function [data_avg,t_avg] = osl_movavg(data,t,winsize,overlap,use_hanning)

% Moving window averaging
% [data_avg,t_avg] = ABMOVAVG(data,t,winsize,overlap)
% -----------------------------------------------------------------
% data        - data to average
% t           - vector of sample times
% winsize     - window to average over (samples)
% overlap     - amount of segment overlapping (e.g. 0.75 = 75% overlap).
%               Default = 0.75.
% use_hanning - apply hanning window to data segment. Default = 0;
% -----------------------------------------------------------------
% AB 2011

if nargin<5;      use_hanning = 0;           end
if use_hanning;   winsize     = winsize *2;  end    
if nargin < 4;    overlap     = 0.75;        end

if use_hanning && overlap > 0,
    warning(['Using a hanning window with overlap may create large edge effects. \n', ...
             '   Test before use. \n']);
end%if

data = [nan(1,winsize),data,nan(1,winsize)];
t    = [nan(1,winsize),t,nan(1,winsize)];

bf       = buffer(data,winsize,round(winsize*overlap),'nodelay');
tbf      = buffer(t,winsize,round(winsize*overlap),'nodelay');

% reject data sections where more than 25% of data is NaN
nans = sum(isnan(bf)) > 0.25*winsize;
bf( :,nans) = nan; % - do this by replacing all with nans
tbf(:,nans) = nan; % - do this by replacing all with nans

if(use_hanning)
hanning_window = repmat(hanning(winsize),1,size(bf,2));
bf = bf.*hanning_window;
end

data_avg = nanmean(bf);
t_avg    = nanmean(tbf);   

data_avg(t_avg==0) = [];
   t_avg(t_avg==0) = [];
   
 data_avg(isnan(t_avg)) = [];
    t_avg(isnan(t_avg)) = [];  
   
end
  
