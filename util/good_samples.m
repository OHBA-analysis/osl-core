function res = good_samples(D,chanind, sampind, trialind)
	% This is a convenience function equivalent to
	%
	% ~any(D.badsamples(chan_inds, sampind, trialind))
	%
	% However, as a modified copy of badsamples(), it is faster
	% because it aggregates over all channels
	%
	% Note - bad channels are always excluded from this check
	% because otherwise the entire timeseries would be bad

	if nargin < 4 || isempty(trialind) 
		trialind = ':';
	end

	if nargin < 3 || isempty(sampind) 
		sampind = ':';
	end

	if nargin < 2 || isempty(chanind) 
		chanind = indchantype(D,'ALL','GOOD'); % By default, check all good channels
	end

	if ischar(sampind) && isequal(sampind, ':')
	    sampind = 1:nsamples(D);
	end

	if ischar(trialind) && isequal(trialind, ':')
	    trialind = 1:ntrials(D);
	end
	

	res = ~badsamples(D, chanind, sampind, trialind);

	% ~any(D.badsamples(chan_inds,:,:))


function res = badsamples(this, chanind, sampind, trialind)
	% Returns an array of 0/1 marking bad data based on artefact events and bad flags
	% FORMAT res = badsamples(this, chanind, sampind, trialind)
	% _______________________________________________________________________
	% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

	% Vladimir Litvak
	% $Id: badsamples.m 6311 2015-01-21 15:44:23Z vladimir $
	res = false(1, nsamples(this), length(trialind));
	chantypes = unique(this.chantype(chanind)); % Channel types to check

	for i = 1:length(trialind)

		% Retrieve all events associated with this trial
	    ev = events(this, trialind(i), 'samples');
	    if iscell(ev)
	        ev = ev{1};
	    end

	    ev = ev(cellfun(@(x) strmatch('artefact',x),{ev.type}) & ismember({ev.value},chantypes)); % These are all the artefact events that apply to the channel types we are inspecting

	    if ~isempty(ev)
	        for k = 1:numel(ev)
	            res(1, ev(k).sample+(0:(ev(k).duration-1)), i) = true;
	        end
	    end
	end

	res = res(:, sampind, :);
	res(:, :, badtrials(this, trialind))  = true;



