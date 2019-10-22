function [metrics,tc] = compute_metrics(D,mains_frequency,artefact_channels)
    % Take in D object, return artefact-related metrics for AFRICA
    %
    % INPUTS
    % - D : MEEG object
    % - mains_frequency : (Hz) mains power will be computed +-1Hz of this freq.
    % - artefact_channels : Compute IC correlations with channels whose
    %                       chantype is in this list
    %
    % OUTPUTS
    % - metrics : struct with various metrics
    % - tc : IC timecourses with NaNs where there were badsamples

    tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';
    num_ics = D.ica.params.num_ics;

    %% COMPUTE METRICS 
    metrics = struct;
    
    % Mains frequency
    tc_fft = fft(demean(tc(:,D.ica.good_samples),2),[],2);
    abs_ft  = abs(tc_fft);
    freq_ax = 0:(D.fsample/2)/(floor(length(abs_ft)/2)-1):(D.fsample/2);
    minhz = 0.1;
    minhz_ind = min(find(freq_ax>minhz));
    spec = abs_ft(1:num_ics,1:floor(length(abs_ft)/2));
    spec(:,1:minhz_ind) = 0;
    [~, ind] = max(spec');
    metrics.mains.mains_frequency = mains_frequency; % Record which mains frequency was used
    metrics.mains.fmax = freq_ax(ind); % Peak frequency for each component
    inds = freq_ax > mains_frequency - 1 & freq_ax < mains_frequency + 1; % Mains +- 1Hz 
    metrics.mains.spec = max(spec(:,inds),[],2); % Maximum power within 1Hz of mains frequency
    spec_n = normalise(spec,2);
    metrics.mains.spec_n = max(spec_n(:,inds),[],2); % Maximum normalised power within 1Hz of mains frequency

    % Kurtosis 
    kurt = kurtosis(tc(:,D.ica.good_samples),[],2);
    kurt = kurt(:);
    metrics.kurtosis.raw = kurt-3;
    metrics.kurtosis.log = log(kurt-3); 
    metrics.kurtosis.abs = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis

    % Cardiac autocorrelation
    bpm_range = [50 100]./60;
    maxlags = 2*D.fsample;
    [~,lags] = xcorr(tc(1,D.ica.good_samples),maxlags,'coeff');
    lags = lags./D.fsample;
    lags_range = lags>bpm_range(1) & lags<bpm_range(2);
    ac_max = zeros(1,num_ics);
    for ic = 1:num_ics
        tc_bp = bandpass(tc(ic,D.ica.good_samples),[0 48],D.fsample);
        ac = xcorr(tc_bp,maxlags,'coeff');
        ac_max(ic) = max(ac(lags_range));
    end
    metrics.cardiac_autocorrelation.value = ac_max(:);

    % Convert chantypes to a list of channel indices to compare to
    chans_to_check = [];
    for j = 1:length(artefact_channels)
        if ischar(artefact_channels{j}) || isstr(artefact_channels{j})
            if isempty(D.indchantype(artefact_channels{j}))
                fprintf(2,'Requested to check correlations with chantype %s, but no channels have this type in the data\n',artefact_channels{j});
            else
                chans_to_check = [chans_to_check D.indchantype(artefact_channels{j})];
            end
        else
            chans_to_check = [chans_to_check artefact_channels{j}];
        end
    end

    tc_bp = osl_filter(tc,[0.1 48],'fs',D.fsample); % filtered ICs
    artefact_data = osl_filter(D(chans_to_check,:),[0.1 48],'fs',D.fsample);

    ac = zeros(size(tc_bp,1),size(artefact_data,1));
    for j = 1:size(tc_bp,1)
        for k = 1:length(chans_to_check)
            artefact_good_samples = ~event_to_sample(D,'artefact_OSL',D.chantype(chans_to_check(1)));
            flt = D.ica.good_samples & artefact_good_samples; % Only use times when both the ICA and the artefact channel were good
            ac(j,k) = abs(nanmedian(osl_movcorr(tc_bp(j,flt)',artefact_data(k,flt)',D.fsample*10,0)));
        end
    end

    for j = 1:length(chans_to_check)
        metric_name = sprintf('corr_chan_%d_%s',chans_to_check(j),D.chantype{chans_to_check(j)});
        metrics.(metric_name).('value')  = ac(:,j);
        metrics.(metric_name).chanind = chans_to_check(j);
    end

    % Assign variables for use in manual GUI
    metrics.mains.value = metrics.mains.spec_n;
    metrics.kurtosis.value = metrics.kurtosis.abs;

    % Add variance as default metric
    metrics.variance.value = nanvar(tc,[],2);
    metrics = orderfields(metrics);
    tc(:,~D.ica.good_samples) = NaN;


