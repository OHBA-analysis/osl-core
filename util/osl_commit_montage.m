function D2 = osl_commit_montage(D,new_filename)
    % Take the current online montage - save an MEEG file with that montage as the 

    if D.montage('getindex') == 0
        fprintf(2,'No montage active - just copying the file\n');
        D2 = D.copy(getfullpath(new_filename));
        return
    else
        m = D.montage('getmontage');
    end


    D2=clone(D.montage('switch',0),getfullpath(new_filename),size(D));

    % Drop all online montages
    D2 = D2.montage('remove',1:D2.montage('getnumber'));

    % Assign the signals;
    D2(:,:,:) = D(:,:,:); 
    D2 = D2.chantype(1:D2.nchannels,D.chantype);
    D2 = D2.chanlabels(1:D2.nchannels,D.chanlabels);
    D2 = D2.units(1:D2.nchannels,D.units);
    D2 = D2.badchannels(1:D2.nchannels,false); % Remove all bad channels
    D2 = D2.badchannels(D.badchannels,true); % Set the bad channels from the input MEEG object

    
    % This is slightly tricky, because we need to re-create the artefacts for the new
    % channel modalities
    D = D.montage('switch',0);
    original_chantypes = D.chantype;
    new_chantypes = D2.chantype;
    new_modalities = unique(new_chantypes);

    for j = 1:D2.ntrials 

        ev = D.events(j);
        
        if iscell(ev)
            assert(length(ev)==1,'Unexpected event length?')
            ev = ev{1};
        end
        
        if ~isempty(ev)
            is_artefact_event = strncmp({ev.type},'artefact',8); % First pick out artefact types - match all possible artefacts
            new_events = ev(~is_artefact_event); % Retain old non-artefacts
            ev = ev(is_artefact_event); % Old artefacts
        end

        % If no artefact events, no need to continue
        if isempty(ev)
            continue
        end

        % Map old modalities to new modalities
        for k = 1:length(new_modalities) % For each new modality
            % Which channels contribute to it?
            tmp = m.tra(strcmp(new_modalities{k},new_chantypes),:); % These are all the rows of tra for the new modality
            old_channel_present = any(tmp,1); % This is whether the old channel contributes to the new modality
            old_modalities = unique(original_chantypes(old_channel_present)); % These are the old modalities
            this_ev = ev(ismember({ev.value},[old_modalities new_modalities{k} 'all'])); % Get all events that match all channels, the new channel, or one of the old channels
            for l = 1:length(this_ev)
                this_ev(l).value = new_modalities{k};
            end
            new_events = [new_events(:);this_ev(:)];
        end
        
        D2 = events(D2,j,new_events);
    end

    D2.save()



