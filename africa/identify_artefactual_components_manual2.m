function [bad_components, metrics, fig_handles,fig_names,fig_titles] = identify_artefactual_components_manual(D,varargin)
    % AfRICA - ArteFact Rejection using Independent Component Analysis
    %
    % Syntax: [bad_components, fig_handles, fig_names, fig_titles] = identify_artefactual_components_manual(S)
    % S needs to contain be produced inside OSL_AFRICA.m
    %
    % OPTIONAL inputs to compute various metrics:
    %   -  ident.do_mains:       compute metric based on mains frequency [0 or 1]
    %   -  ident.do_kurt:        compute metric based on extreme kurtosis [0 or 1]
    %   -  ident.do_cardiac:     compute metric based on similarity to cardiac cycle [0 or 1]
    %   -  ident.artefact_channels: compute metric based on correlation with external channels
    %                               can be either part of a channame or a specific
    %                               channel idx [e.g. {'ECG','EOG','312'}]
    %   -  ident.launch_gui: if false, compute & save metrics & topos but don't use GUI to choose bad components
    % Output:
    %   - bad_components: A list of the components identified as bad.
    %
    % Romesh Abeysuriya 2017
    % Adam Baker 2014

    arg = inputParser;
    arg.addParameter('do_mains',false); 
    arg.addParameter('mains_frequency',50); 
    arg.addParameter('do_kurt',false); 
    arg.addParameter('do_cardiac',false); 
    arg.addParameter('artefact_channels',[]); %  cellstr(args.artefact_channels);? 
    arg.addParameter('launch_gui',true); 
    arg.parse(varargin{:});
    args = arg.Results;

    if args.do_cardiac
        if any(strcmpi(args.artefact_channels,'ECG'))
            do_cardiac_autocorrelation = 0;
        else
            do_cardiac_autocorrelation = 1;
        end
    else
        do_cardiac_autocorrelation = 0;
    end

    fig_handles = [];
    fig_names   = [];
    fig_titles  = [];
    bad_components = []; % If component identification isn't used, don't return any bad components
    %% Select only good data for classification

    % Good samples/trials
    samples_of_interest = ~all(badsamples(D,':',':',':'));
    samples_of_interest = reshape(samples_of_interest,1,D.nsamples*D.ntrials);

    sm = D.ica.sm;
    tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';

    num_ics = D.ica.params.num_ics;

    tc_fft = fft(demean(tc(:,samples_of_interest),2),[],2);

    abs_ft  = abs(tc_fft);
    freq_ax = 0:(D.fsample/2)/(floor(length(abs_ft)/2)-1):(D.fsample/2);

    minhz = 0.1;
    minhz_ind = min(find(freq_ax>minhz));

    spec = abs_ft(1:num_ics,1:floor(length(abs_ft)/2));
    spec(:,1:minhz_ind) = 0;

    %% COMPUTE METRICS 
    metrics = struct;
    
    % Mains frequency
    if args.do_mains
        spec_n = normalise(spec,2);
        inds = freq_ax > args.mains_frequency - 1 & freq_ax < args.mains_frequency + 1;
        metrics.mains.value = max(spec_n(:,inds),[],2);
    end

    % Kurtosis 
    if args.do_kurt
        kurt = kurtosis(tc(:,samples_of_interest),[],2);
        kurt_t = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis
        metrics.kurtosis.value = kurt_t(:);
    end

    % Cardiac  
    if do_cardiac_autocorrelation
        disp('No ECG artefact channels specified - detecting cardiac components without ECG')

        bpm_range = [50 100]./60;
        maxlags = 2*D.fsample;

        [~,lags] = xcorr(tc(1,samples_of_interest),maxlags,'coeff');
        lags = lags./D.fsample;
        lags_range = lags>bpm_range(1) & lags<bpm_range(2);

        ac_max = zeros(1,num_ics);
        for ic = 1:num_ics
            tc_bp = bandpass(tc(ic,samples_of_interest),[0 48],D.fsample);
            ac = xcorr(tc_bp,maxlags,'coeff');
            ac_max(ic) = max(ac(lags_range));
        end

        metrics.cardiac_autocorrelation.value = ac_max(:);
    end

    % Detect artefact channels related components
    if ~isempty(args.artefact_channels)

        for artefact_chantype = unique(args.artefact_channels)
            if ~isfield(metrics,artefact_chantype)

                if isempty(str2num(cell2mat(artefact_chantype)))
                    % We have a string, extract relevant channels matching channel type name
                    artefact_data = D(find(strcmp(D.chantype,artefact_chantype)),:);
                    artefact_data = reshape(artefact_data,size(artefact_data,1),[]);
                    metric_name = cell2mat(artefact_chantype);
                else
                    % We have a channel number
                    artefact_data = D(str2num(cell2mat(artefact_chantype)),:);
                    metric_name = ['Chan_' cell2mat(artefact_chantype)];
                end

                % Bandpass filter
                for ac = 1:size(artefact_data,1)
                    artefact_data(ac,:) = bandpass(artefact_data(ac,:),[0.1 48],D.fsample);
                end

                ac_corr = zeros(1,num_ics);
                for ic = 1:num_ics
                    tc_bp = bandpass(tc(ic,samples_of_interest),[0.1 48],D.fsample);
                    %ac_corr(ic) = sum(abs(corr(tc_bp',artefact_data')));
                    for ac = 1:size(artefact_data,1)
                        ac_corr(ic) = ac_corr(ic) + abs(nanmedian(osl_movcorr(tc_bp',artefact_data(ac,find(samples_of_interest))',D.fsample*10,0)));
                    end
                end
                ac_corr = ac_corr ./size(artefact_data,1);
                
                % output correlation and timecourse
                metrics.(metric_name).('value')  = ac_corr(:);
                metrics.(metric_name).timeCourse = artefact_data.';
                metrics.(metric_name).timeCourse(~samples_of_interest,:) = nan;
            end
        end

    end

    % Add variance as default metric;
    metrics.variance.value = var(sm);
    metrics = orderfields(metrics,[find(strcmp(fieldnames(metrics),'variance'));find(~strcmp(fieldnames(metrics),'variance'))]);
    tc(:,~samples_of_interest) = nan;

    if args.launch_gui
        bad_components = identify_artefactual_components_manual_gui(D,tc,D.ica.topos,metrics,D.ica.bad_components);
    end
