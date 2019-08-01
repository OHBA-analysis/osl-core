function [out,cid,sid,tid] = good_samples(D, cid, sid, tid, check_ctype)
%
% [out,chanind,sampind,trialind] = good_samples(D, chanind, sampind, trialind, check_chantype)
%
% This is a convenience function equivalent to
%
%   ~any(D.badsamples(chanind, sampind, trialind),1)
%
% i.e. the outputs should match exactly. Behaviours do differ
% slightly though if bad channels are included in 
% chanind because they are ignored by D function. Similarly,
% bad channel types are excluded when finding good samples. Basically
% D function should do what you expect if you are not using the
% data from bad channels.
%
% However, as a modified copy of badsamples(), it is faster
% because it aggregates over all channels
%
% Romesh Abeysuriya 2017, JH 2019

    if nargin < 5 
        check_ctype = true; 
    end
    
    if nargin < 4 || isempty(tid) 
		tid = ':';
    end

    if nargin < 3 || isempty(sid) 
		sid = ':';
    end

    if nargin < 2 || isempty(cid) 
		cid = D.indchantype('ALL','GOOD'); % By default, check all good channels
	elseif ischar(cid) && isequal(cid, ':')
		cid = D.indchantype('ALL');
    elseif ischar(cid)
        cid = D.indchantype(cid);
    elseif iscellstr(cid)
        cid = D.indchantype(cid{:});
    elseif islogical(cid)
        assert( numel(cid)==D.nchannels, 'Bad channel-mask size.' );
        cid = find(cid);
    end

    if ischar(sid) && isequal(sid, ':')
	    sid = 1:nsamples(D);
    elseif islogical(sid)
        assert( numel(sid)==D.nsamples, 'Bad sample-mask size.' );
        sid = find(sid);
    end

    if ischar(tid) && isequal(tid, ':')
	    tid = 1:ntrials(D);
    elseif ischar(tid)
        tid = D.indtrial(tid);
    elseif iscellstr(tid)
        tid = D.indtrial(tid{:});
    elseif islogical(tid)
        assert( numel(tid)==D.ntrials, 'Bad trial-mask size.' );
        tid = find(tid);
    end

    if isempty(intersect( D.badchannels, cid ))
        out = true(1, nsamples(D), length(tid));
    else
        out = false(1, numel(sid), length(tid));
        return;
    end
	chantypes = unique(D.chantype(cid)); % Channel types to check

	% If online montage, we also need to check the constituent channels
    if D.montage('getindex') > 0 
		m = D.montage('getmontage');
		Dtemp = D.montage('switch',0);
        old_chanind = find(any(m.tra(cid,:),1));
        old_chantypes = unique(Dtemp.chantype(old_chanind)); %#ok
        chantypes = union(chantypes,old_chantypes);
    end
    
    chantypes{end+1} = 'all'; % Make sure we support modality-independent artefacts
    chantypes = upper(chantypes);
    
	for i = 1:length(tid)

		% Retrieve all events associated with D trial
                
	    ev = events(D, tid(i), 'samples');	            
        if iscell(ev), ev = ev{1}; end
        if isempty(ev), continue; end
        
	    ev = ev(strncmp( {ev.type}, 'artefact', 8 )); % First pick out artefact types
	    if isempty(ev), continue; end
	    
        % Then filter by chantype. 
        % This prevents bugs if the trial event value is not a chantype but a number
        if check_ctype 
            ev = ev(ismember(upper({ev.value}),chantypes)); 
        end
    
	    if ~isempty(ev)
	        for k = 1:numel(ev)
	            out(1, ev(k).sample+(0:(ev(k).duration-1)), i) = false;
	        end
	    end

	end

	out = out(:, sid, :);
	out(:, :, badtrials(D, tid)) = false;

end