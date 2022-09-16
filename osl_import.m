function D = osl_import(raw_data,varargin)
    % [D fname] = osl_import(raw_data,varargin)
    % 
    % Converts .fif file to continuous spm data file
    %
    % INPUTS
    % - raw_data : file name of input .fif file or .ds folder
    % - varargin : see below
    %
    % OUTPUTS
    % - D : MEEG object
    %
    % RA 2017
    % MW

    arg = inputParser;
    arg.addParameter('outfile',[],@ischar); % Output file name
    arg.addParameter('trigger_channel_mask','0000000001111111'); % !?!
    arg.addParameter('artefact_channels',[]); % Specify artefact channels that will be retained based on chantype
    arg.addParameter('other_channels',[]); % Specify other channels to retain based on chanlabel
    arg.addParameter('num_condition_trials_expected',[]); % number of expected trials for each condition
    arg.addParameter('condition_code',[]); % list of condition trigger codes 
    arg.addParameter('bad_segments',[]); % n x 2 matrix of bad times to add to the output MEEG object 
    arg.addParameter('bad_segment_modalities',{}); % Cell array of modalities to apply bad segments to - default is all MEGANY modalities that are present

    arg.parse(varargin{:});
    S = arg.Results;

    [pathstr,fname,ext] = fileparts(raw_data);
    is_ctf = strcmp(ext,'.ds');

    if isempty(S.outfile)
        S.outfile = fullfile(pathstr,[fname '.mat']);
    end

    % Create directory if it doesn't exist
    output_path = fileparts(S.outfile);
    if ~isdir(output_path)
        mkdir(output_path);
    end

    if is_ctf
        %S.channels=ctf_chans;     % HL Mod 1.2 - CTF channel selection
        %TC 2013 the above works only for a 275 channel ctf system
        %To make more general, we want to read stim channels, the clock channel
        %and meg grads. Ideally data should actually be beamformed in 3rd gradient
        %which would require the reference channels but that would require
        %modification of the forward model (and you would also need to read in
        %balancing coefs). For now the data is read in unbalanced (raw) - this
        %does not take advantage of the noise corrections available to CTF
        %systems.
        S.channels='ctf_grads';    
        S.checkboundary = 0;  
    else
        S.channels = 'all';
        S.checkboundary = 1;
        S.fixoxfordneuromag = S.trigger_channel_mask;
    end

    S.dataset = raw_data;
    S.usetrials = 1;
    S.datatype = 'float32-le';
    S.eventpadding = 0;
    S.saveorigheader = 1;
    S.conditionlabel = {'Undefined'};
    S.inputformat = [];
    S.mode='continuous';

    D = spm_eeg_convert_4osl(S);

    %% Add coil to channel mapping in D.sensors (used later for local spheres)
    if ~isempty(D.sensors('MEG'))
        sens = D.sensors('MEG');
        sens.coilchan = sens.tra==1;
        D = sensors(D,'MEG',sens);
    end

    %% If any events have value that is not > 0, set it to 1
    %% MWW added back in from osl1.2.beta.15 - Oct 2013
    ev = events(D,1);
    for j = 1:length(ev)
        try
            if isempty(ev(j).value) || ~(ev(j).value>0)
                ev(j).value=1;
            end
        catch
            warning('Strange event type encountered in the raw data object. Please inspect the events in the raw file!');
        end
    end
    D = D.events(1,ev);

    % Check to see if correct number of events have been found
    if ~is_ctf
        if ~isempty(ev) 
            vals = [ev.value];
            fprintf('Found %d events\n',length(vals));

            if ~isempty(S.num_condition_trials_expected)
                for j = 1:length(S.condition_code)
                    n_events = sum(vals==S.condition_code(j));
                    if n_events ~= S.num_condition_trials_expected(j)
                        fprintf(2,'Warning - expected %d events with value %d, but found %d\n',S.num_condition_trials_expected(j),S.condition_code(j),n_events);
                    end
                end
            end

        else
            fprintf(2,'No events detected\n');
        end
    end

    if ~isempty(S.bad_segments)
        if isempty(S.bad_segment_modalities)
            S.bad_segment_modalities = unique(D.chantype(D.indchantype('MEGANY')));
        end
        BadEvents = [];
        for j = 1:length(S.bad_segment_modalities)
            for k = 1:size(S.bad_segments,1)
                BadEvents(end+1).type   = 'artefact_OSL';
                BadEvents(end).value  = S.bad_segment_modalities{j};
                BadEvents(end).time   = S.bad_segments(k,1);
                BadEvents(end).duration = diff(S.bad_segments(k,:));
                BadEvents(end).offset = 0;
            end
        end
        ev = D.events;
        D = D.events(1,[ev(:); BadEvents(:)]);
    end

    D.save;

    chantypes = {'EOG','ECG'};
    for j = 1:length(chantypes)
        if isempty(D.indchantype(chantypes{j}))
            fprintf(2,'No %s channels were automatically identified - you may need to specify them manually with D=D.chantype(XXX,''%s'') if you wish to use them e.g. in AFRICA\n',chantypes{j},chantypes{j})
        end
    end
