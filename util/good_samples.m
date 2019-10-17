function res = good_samples(D,chanind, sampind, trialind, checkchantype)
    % This is a convenience function equivalent to
    %
    % ~any(D.badsamples(chanind, sampind, trialind))
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
    % Romesh Abeysuriya 2017

    if nargin < 5 
        checkchantype = 1; 
    end
    
    if nargin < 4 || isempty(trialind) 
        trialind = ':';
    end

    if nargin < 3 || isempty(sampind) 
        sampind = ':';
    end

    if nargin < 2 || isempty(chanind) 
        chanind = D.indchantype('ALL','GOOD'); % By default, check all good channels
    elseif ischar(chanind) && isequal(chanind, ':')
        chanind = D.indchantype('ALL');
    end

    if ischar(sampind) && isequal(sampind, ':')
        sampind = 1:nsamples(D);
    end

    if ischar(trialind) && isequal(trialind, ':')
        trialind = 1:ntrials(D);
    end
    

    res = true(1, nsamples(D), length(trialind));
    chantypes = unique(D.chantype(chanind)); % Channel types to check

    % If online montage, we also need to check the constituent channels
    if D.montage('getindex') > 0 
        m = D.montage('getmontage');
        Dtemp = D.montage('switch',0);
        old_chantypes = unique(Dtemp.chantype(find(any(m.tra(chanind,:),1))));
        if ~isempty(old_chantypes)
            chantypes = union(chantypes,old_chantypes);
        end
    end
    
    chantypes{end+1} = 'all'; % Make sure we support modality-independent artefacts

    for i = 1:length(trialind)

        % Retrieve all events associated with D trial
                
        ev = events(D, trialind(i), 'samples');             
        if iscell(ev)
            ev = ev{1};
        end
        
        if isempty(ev)
            continue
        end
        
        ev = ev(strncmp({ev.type},'artefact',8)); % First pick out artefact types

        if isempty(ev)
            continue
        end
        
        if checkchantype % Then filter by chantype. This prevents bugs if the trial event value is not a chantype but a number
            ev = ev(ismember(upper({ev.value}),upper(chantypes))); 
        end
    
        if ~isempty(ev)
            for k = 1:numel(ev)
                res(1, ev(k).sample+(0:(ev(k).duration-1)), i) = false;
            end
        end

    end

    res = res(:, sampind, :);
    res(:, :, badtrials(D, trialind))  = false;

    bad_channels = intersect(D.badchannels,chanind);
    if any(bad_channels)
        res(:) = false;
    end


