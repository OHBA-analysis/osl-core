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
    arg.addParameter('modalities',[]); 
    arg.addParameter('max_iter',3);
    arg.addParameter('max_bad_channels',10); % Maximum number of bad channels to identify at this stage 
    arg.addParameter('badchannels',true); % Check for bad channels
    arg.addParameter('badtimes',true); % Check for bad events
    arg.addParameter('dummy_epoch_tsize',1); % Dummy epoch size in seconds
    arg.addParameter('measure_fns',{'std'}); % list of outlier metric func names to use for bad segment marking
    arg.addParameter('event_threshold',0.3); % list of robust GLM weights thresholds to use on EVs for bad segment marking, the LOWER the theshold the less aggressive the rejection
    arg.addParameter('channel_threshold',0.01); % list of robust GLM weights thresholds to use on chans for bad segment marking, the LOWER the theshold the less aggressive the rejection
    arg.parse(varargin{:});
    options = arg.Results;

    D_original = D;
    
    if isempty(options.modalities)
        options.modalities = unique(D.chantype(D.indchantype({'MEG','MEGGRAD','MEGPLANAR','EEG'})))
        fprintf('Detecting artefacts in channel types: %s\n',strjoin(options.modalities,','));
    end
        
    % Are we continuous and checking for bad timess? If continuous, D becomes
    % a temporary copy with artificial epochs. Otherwise, D is just the
    % original input file
    continuous = strcmp(D_original.type,'continuous');
    if continuous && arg.Results.badtimes

        tind_tsize = arg.Results.dummy_epoch_tsize*D_original.fsample;
        start_tinds=1:tind_tsize:(D_original.nsamples);
        Ndummytrials=length(start_tinds);

        epochinfo.conditionlabels=cell(Ndummytrials-1,1);
        epochinfo.trl=zeros(Ndummytrials-1, 3);

        for ii=1:Ndummytrials-1,
            epochinfo.conditionlabels{ii}='Dummy';
            epochinfo.trl(ii,1)=start_tinds(ii); % trl is in time indices
            epochinfo.trl(ii,2)=min(start_tinds(ii)+tind_tsize-1);
            epochinfo.trl(ii,3)=0;
        end

        S = epochinfo;
        S.D = D_original;
        D = osl_epoch(S);
        D = D.montage('switch',D_original.montage('getindex'));
    else
        D = D_original;
    end
    
    
    % For each modality, detect badness
    for mm = 1:length(arg.Results.modalities)

        modality=arg.Results.modalities{mm};
        
        chan_list = find(strcmp(D.chantype(),modality));
        chanind = setdiff(chan_list, D.badchannels); % All currently good chans

        iters = 0;
        detected_badness = true; % The while loop continues as long as something bad was found
        
        while detected_badness && iters < arg.Results.max_iter
            iters = iters+1;
            detected_badness = false; % Terminate by default unless something bad is found

            if arg.Results.badchannels % Do a pass for bad channels
                for ii = 1:length(arg.Results.measure_fns)

                    dat = D(chanind,:,D.indtrial(D.condlist,'good')); % Only work on trials that haven't been rejected at any point
                    datchan = feval(arg.Results.measure_fns{ii},reshape(dat,size(dat,1),size(dat,2)*size(dat,3)),[],2);
                    
                    [b,stats] = robustfit(ones(length(datchan),1),datchan,'bisquare',4.685,'off');

                    if numel(arg.Results.channel_threshold) == 1
                        channel_threshold = arg.Results.channel_threshold;
                    else
                        channel_threshold = arg.Results.channel_threshold(ii);
                    end

                    [bad_chan] = find(stats.w < channel_threshold);
                    bad_chan_w = stats.w(bad_chan);
                    [sorted_w iw] = sort(bad_chan_w);
                    sorted_bad_chan = bad_chan(iw);

                    if length(badchannels(D))+length(sorted_bad_chan) > arg.Results.max_bad_channels
                        warning(['More than arg.Results.max_bad_channels=' num2str(arg.Results.max_bad_channels) ' have been detected. But only marking the worst ' num2str(arg.Results.max_bad_channels) ' as bad']);
                    end

                    num_to_add=min(length(sorted_bad_chan),arg.Results.max_bad_channels-length(badchannels(D)));
                    
                    if num_to_add > 0
                        sorted_bad_chan = sorted_bad_chan(1:num_to_add);
                    else
                        sorted_bad_chan = [];
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

            if arg.Results.badtimes % Do a pass for bad trials
                for ii = 1:length(arg.Results.measure_fns)
                    trials = D.indtrial(D.condlist,'good');
                    dat = D(chanind,:,trials);
                    datchan = feval(arg.Results.measure_fns{ii},reshape(dat,size(dat,1)*size(dat,2),size(dat,3)),[],1);
                    [b,stats] = robustfit(ones(length(trials),1),datchan,'bisquare',4.685,'off');

                    if numel(arg.Results.event_threshold) == 1
                        event_threshold=arg.Results.event_threshold;
                    else
                        event_threshold=arg.Results.event_threshold(ii);
                    end

                    bad_ev = find(stats.w < event_threshold);

                    if length(bad_ev) > 0
                        rej = zeros(1,D.ntrials);
                        rej(D.badtrials) = 1;
                        rej(trials(bad_ev)) = 1;
                        D = D.badtrials(1:length(rej), rej);  
                        detected_badness = true;
                    end

                    % Mark bad trials for exclusion at the next iteration
                    trials = D.indtrial(D.condlist,'good');
                end
            end
        end
    end

    % If we made a temporary epoched file, need to write the bad
    % segments back to the original file for output
    if continuous && arg.Results.badtimes
        % Initialize output MEEG
        D_out = D_original;

        % Copy bad channels from temporary D 
        bad_channels = D.badchannels;
        if ~isempty(bad_channels)
            D_out = D_out.badchannels(bad_channels,1);
        end

        % Copy bad segments from temporary D
        BadEpochs = {};
        badtrialstmp = D.badtrials;
        badtrials = zeros(1,D.ntrials);
        badtrials(badtrialstmp)=1;

        diffbadtrials = diff([0 badtrials]);
        ups=find(diffbadtrials == 1);
        downs=find(diffbadtrials == -1);
        if ~isempty(ups)

            if(badtrials(end))
                downs(end+1) = length(badtrials);
            end

            for jj = 1:length(ups)
                BadEpochs{jj}(1) = D_out.time(epochinfo.trl(ups(jj),1));
                BadEpochs{jj}(2) = D_out.time(epochinfo.trl(downs(jj),2));
            end

            D_out = set_bad(D_out,BadEpochs);

        end

        % Remove the temporary epoched D
        % DANGER - Fix this - command below is only safe if the if condition in this 
        % block is exactly the same as the one that makes the dummy MEEG at the start
        % of the file!
        D.delete();
    else
        D_out = D;
    end

    bc = D_out.badchannels;
    for j = 1:length(bc)
        fprintf('Channel %d (%s) is bad\n',bc(j),D.chanlabels{bc(j)});
    end

    if continuous
        ev = D_out.events;
        for j = 1:numel(ev)
            if isfield(ev,'type') && strcmp(ev(j).type,'artefact_OSL')
                fprintf('Bad epoch from %.2f-%.2f\n',ev(j).time,ev(j).time + ev(j).duration)
            end
        end
    else
        bt = D_out.badtrials;
        for j = 1:length(bt)
            fprintf('Trial %d is bad\n',bt(j));
        end
    end

    D = D_out;
    
function D = set_bad(D,BadEpochs)
    % Save bad epochs using method meeg/events
    BadEvents = struct([]);
    for ev = 1:numel(BadEpochs)
        if numel(BadEpochs{ev} == 2)
            BadEvents(ev).type     = 'artefact_OSL';
            BadEvents(ev).value    = 'all';
            BadEvents(ev).time     =  BadEpochs{ev}(1);
            BadEvents(ev).duration = diff(BadEpochs{ev});
            BadEvents(ev).offset = 0;
        end
    end

    % Load events
    Events = D.events;

    % Remove previous bad epoch events
    if isfield(Events,'type')
        Events(strcmp({Events.type},'artefact_OSL')) = [];
    end

    % Concatenate new and old events
    if size(Events,1) < size(Events,2)
        BadEvents = BadEvents(:);
    end

    if ~isempty(BadEvents)
      Events = [Events(:); BadEvents(:)];
    end

    % Merge new events with previous
    D = events(D,1,Events);
