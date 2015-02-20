function c = osl_movcorr(s1,s2,winsize,overlap,resamp)
% Moving window correlation
% c = OSL_MOVCORR(s1,s2,winsize,overlap,resamp)
% -----------------------------------------------------------------
% s1,s2   - [1 x Nsamples] signals to calculate correlation between
% winsize - window to average correlation over [samples]
% overlap - amount of segment overlapping (e.g. 0.75 = 75% overlap)
% -----------------------------------------------------------------
% AB 2011

if ~exist('overlap','var') || isempty(overlap)
  overlap = 0.75;
end

if ~exist('resamp','var') || isempty(resamp) 
  resamp = 0;
end

if resamp == 1;
  resamp = 'nearest';
elseif resamp == 0;
  resamp = [];
end
  
% Chunk data into overlapping sections
s1w = buffer(s1,winsize,round(winsize*overlap),'nodelay');
s2w = buffer(s2,winsize,round(winsize*overlap),'nodelay');

s1w(s1w(:,end)==0,end) = nan;
s2w(s2w(:,end)==0,end) = nan;

% Compute correlation between corresponding sections
c = zeros(size(s1w,2),1);
for i = 1:size(s1w,2)
  % Find NaNs
  notNaN =  ~any(isnan([s1w(:,i),s2w(:,i)]),2);
  % If NaNs make up more than 25% of the window then return NaN, otherwise
  % compute the correlation excluding these values.
  if sum(notNaN) < 0.75*size(s1w,1)
    notNaN(:) = 1; % this ensures NaN values are passed to prcorr which will return NaN.
  end  
  try
      c(i) = prcorr2(s1w(notNaN,i),s2w(notNaN,i));
  catch
      c(i) = corr(s1w(notNaN,i),s2w(notNaN,i));
  end
  %plot(s1w(:,i)); title(num2str(c(i)));
end

if ~isempty(resamp)
  
  % Get start and end limits, accounting for zero padding at end of signal
  n_pad = find(cumsum(s1w(end:-1:1,end))==0, 1, 'last' );
  if isempty(n_pad)
    n_pad = 0;
  end
  t_start = round(winsize/2);
  t_end = length(s1) - round(winsize/2) + n_pad;
  
  
  
  % Resample to same size as original data
  ti = t_start:1:t_end;
  t = linspace(t_start,t_end,length(c));
  ci = nan(size(s1));
  warning('off', 'MATLAB:interp1:NaNinY');
  ci(ti) = interp1(t,c,ti,resamp);
  warning('on', 'MATLAB:interp1:NaNinY');
  
  % Remove excess at end due to padding
  ci(end-round(winsize/2)+1:end) = NaN; 
  
  c = ci;
end

end
 
