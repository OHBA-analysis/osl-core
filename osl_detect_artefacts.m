function D = osl_detect_artefacts(D,varargin)
    % Detect artefacts in MEG sensor data
    %
    % This is an adapted version of osl_detect_badevent which checks for both bad epochs/trials and
    % bad channels. The basic workflow is to detect bad channels, and bad times (which are segments for continuous
    % data, and trials for data that has been epoched). 
    % 
    % The input is an MEEG, and the output is a modified version of the input MEEG, not saved to disk, with
    % the 'badchannels' and either 'badtrials' updated, or a set of 'events' of type 'artefact_OSL' added. 
    % Optional inputs let you set whether to check for bad channels and/or bad times, and what outlier function and
    % thresholds to use. 
    %
    % Bad segments are identified in continuous recordings by segmenting into
    % artificial trials, identifying bad trials, and then mapping those 
    % trials back onto the continuous recording. So the workflow is
    %
    % CONTINUOUS - Badchannels and Badsegments
    %   - Check for bad channels
    %   - Segment into trials
    %   - Check for bad channels and bad trials
    %   - Map bad trials back to segments of continuous recording
    %
    % Whereas for epoched data, we just look for bad trials directly
    %
    % EPOCHED - Badchannels and Badtrials
    %   - Check for bad channels and bad trials
    % 
    % Could consider it as - check for badchannels or badtimes
    %
    % PROGRAM LOGIC
    %
    % Continuous - badchannels=true,badtimes=false
    %   - Detect bad channels in original D object
    % Continuous - badchannels=true/false,badtimes=true
    %   - Split into artificial epochs
    %   - Detect bad events and optionally bad channels
    %   - Map events back to times in continuous recording
    % Epoched - badchannels=true,badtimes=true/false
    %   - Detect bad channels and/or bad events in original D object
    %
    % RA 2017
    % MWW 2013

    arg = inputParser;
    arg.addParameter('modalities',{},@iscell); 
    arg.addParameter('max_iter',3);
    arg.addParameter('max_bad_channels',10); % Maximum number of bad channels allowed after this function returns - including any already marked bad
    arg.addParameter('badchannels',true); % Check for bad channels
    arg.addParameter('badtimes',true); % Check for bad events
    arg.addParameter('dummy_epoch_tsize',1); % Dummy epoch size in seconds
    arg.addParameter('measure_fns',{'std'}); % list of outlier metric func names to use for bad segment marking
    arg.addParameter('event_threshold',0.3); % list of robust GLM weights thresholds to use on EVs for bad segment marking, the LOWER the theshold the less aggressive the rejection
    arg.addParameter('channel_threshold',0.01); % list of robust GLM weights thresholds to use on chans for bad segment marking, the LOWER the theshold the less aggressive the rejection
    arg.addParameter('artefact_type_name',[]); % Specify a name for the artefacts to be added. All existing artefacts with this name will be discarded
    arg.parse(varargin{:});
    options = arg.Results;
    
    if isempty(options.modalities)
        options.modalities = unique(D.chantype(D.indchantype({'MEG','MEGMAG','MEGGRAD','MEGPLANAR','EEG'})));
        fprintf('Detecting artefacts in channel types: %s\n',strjoin(options.modalities,','));
    end

    if isempty(options.artefact_type_name) && any(~ismember(options.modalities,{'MEG','MEGMAG','MEGGRAD','MEGPLANAR','EEG'}))
        fprintf(2,'Modalities contains non-MEG channels, so artefacts are being automatically marked as ''OSL_non_meg_artefact'' so they will not interact with badsamples, but will be used in AFRICA when computing artefact correlations\n');
        fprintf(2,'You can set the ''artefact_type_name'' to override this.\n');
        options.artefact_type_name = 'OSL_non_meg_artefact';
    elseif isempty(options.artefact_type_name)
        options.artefact_type_name = 'artefact_OSL';
    end


    % Are we continuous and checking for bad times? If continuous, D becomes
    % a temporary copy with artificial epochs. Otherwise, D is just the
    % original input file
    continuous = strcmp(D.type,'continuous');

    if continuous
        existing_badsamples = event_to_sample(D,options.artefact_type_name); % Find samples that are already bad *based on the same event type as the new proposed badness*
        dummy_trialsize = options.dummy_epoch_tsize*D.fsample;
        dummy_ntrials = floor(D.nsamples/(options.dummy_epoch_tsize*D.fsample));
        convert_to_trial = @(x) reshape(x(:,1:dummy_trialsize*dummy_ntrials),size(x,1),dummy_trialsize,dummy_ntrials); % This function 'epochs' a matrix from chans x continuous_time to chans x trial_time x trials
        data = convert_to_trial(D(:,:));
        badtrials = squeeze(any(convert_to_trial(existing_badsamples),2)); % Existing bad trials
    else
        data = D(:,:,:);
        badtrials = ismember(1:D.ntrials,D.badtrials); % Currently bad trials
    end
    
    trials_to_reject = struct; % For each modality, store which trials were rejected

    % For each modality, detect badness
    for mm = 1:length(options.modalities)

        modality=options.modalities{mm};
        
        chan_list = find(strcmp(D.chantype(),modality));
        chanind = setdiff(chan_list, D.badchannels); % All currently good chans

        iters = 0;
        detected_badness = true; % The while loop continues as long as something bad was found
        
        while detected_badness && iters < options.max_iter
            iters = iters+1;
            detected_badness = false; % Terminate by default unless something bad is found

            if options.badchannels && length(chanind) > 1  % Do a pass for bad channels as long as >1 channel remains
                for ii = 1:length(options.measure_fns)

                    dat = data(chanind,:,~badtrials); % Only work on trials that haven't been rejected at any point
                    datchan = feval(options.measure_fns{ii},reshape(dat,size(dat,1),size(dat,2)*size(dat,3)),[],2);
                    
                    [b,stats] = robustfit(ones(length(datchan),1),datchan,'bisquare',4.685,'off');

                    if numel(options.channel_threshold) == 1
                        channel_threshold = options.channel_threshold;
                    else
                        channel_threshold = options.channel_threshold(ii);
                    end

                    bad_chan = find(stats.w < channel_threshold);
                    [~,iw] = sort(stats.w(bad_chan));
                    sorted_bad_chan = bad_chan(iw);
                    to_add = setdiff(sorted_bad_chan,D.badchannels); % Remove already marked bad channels
                    max_add = options.max_bad_channels-length(D.badchannels); % Maximum number of channels to add
                    if length(to_add) > max_add
                        fprintf(2,'Detected %d new bad channels, only rejecting the worst %d so there will be options.max_bad_channels=%d total\n', length(to_add), max_add,options.max_bad_channels);
                        sorted_bad_chan = sorted_bad_chan(1:max_add);
                    end
                    
                    % set bad channels in D
                    if length(sorted_bad_chan) > 0              
                        bad_channels=chanind(sorted_bad_chan);
                        D = badchannels(D, bad_channels, ones(length(bad_channels),1));
                        detected_badness = true;
                    end

                    % correct chan inds for next bit
                    chanind = setdiff(chan_list, D.badchannels);
                end
            end

            if options.badtimes % Do a pass for bad trials
                for ii = 1:length(options.measure_fns)
                    goodtrials = find(~badtrials);
                    dat = data(chanind,:,goodtrials);
                    datchan = feval(options.measure_fns{ii},reshape(dat,size(dat,1)*size(dat,2),size(dat,3)),[],1);
                    [b,stats] = robustfit(ones(length(datchan),1),datchan,'bisquare',4.685,'off');

                    if numel(options.event_threshold) == 1
                        event_threshold=options.event_threshold;
                    else
                        event_threshold=options.event_threshold(ii);
                    end

                    % For each good trial, mark that trial number as bad
                    new_badtrials = find(stats.w < event_threshold);
                    badtrials(goodtrials(new_badtrials))=1;
                end
            end
        end
    end

    if options.badtimes

        if continuous

            % Find the onset and offset times
            db = find(diff([0;badtrials]));
            onset = db(1:2:end);
            offset = db(2:2:end);
            if length(offset)<length(onset) % This means that the state changed without an end of change
                offset(end+1) = length(badtrials)+1; % It overruns by 1 at the end
            end

            BadEvents = struct([]);
            for j = 1:length(onset)
                BadEvents(j).type     = options.artefact_type_name;
                BadEvents(j).value    = 'all';
                BadEvents(j).time     =  (onset(j)-1)*options.dummy_epoch_tsize+ 1/D.fsample;
                BadEvents(j).duration =  (offset(j)-onset(j))*options.dummy_epoch_tsize;
                BadEvents(j).offset = 0;
            end

            Events = D.events;

            % Remove previous bad epoch events of this same type - they are already accounted for at the start
            if isfield(Events,'type')
                Events(strcmp({Events.type},options.artefact_type_name)) = [];
            end

            % Concatenate new and old events
            if size(Events,1) < size(Events,2)
                BadEvents = BadEvents(:);
            end

            if ~isempty(BadEvents)
              Events = [Events(:); BadEvents(:)];
            end

            % Merge new events with previous
            D = D.events(1,Events);
        else
            % Map directly to bad trials
            D = D.badtrials(1:length(badtrials),badtrials);
        end
        
    end

    % Display summary at the end
    bc = D.badchannels;
    for j = 1:length(bc)
        fprintf('Channel %d (%s) is bad\n',bc(j),D.chanlabels{bc(j)});
    end

    if continuous
        ev = D.events;
        ev_types = unique({ev.type});
        for j = 1:length(ev_types)
            fprintf('BAD TIMES - EVENT TYPE "%s"\n',ev_types{j});
            for k = 1:length(ev)
                if strcmp(ev(k).type,ev_types{j})
                    fprintf('Bad epoch from %.2f-%.2f\n',ev(k).time,ev(k).time + ev(k).duration)
                end
            end
        end
    else
        bt = D.badtrials;
        for j = 1:length(bt)
            fprintf('Trial %d is bad\n',bt(j));
        end
    end
