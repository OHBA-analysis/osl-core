function [metrics,tc] = compute_metrics(D,mains_frequency,artefact_channels)
	% Take in D object, return artefact-related metrics
	%
	% INPUTS
	% - D : MEEG object
	% - mains_frequency : (Hz) mains power will be computed +-1Hz of this freq.
	% - artefact_channels : Compute IC correlations with channels whose
	%   					chantype is in this list
	%
	% OUTPUTS
	% - metrics : struct with various metrics
	% - tc : IC timecourses with NaNs where there were badsamples

	samples_of_interest = ~all(badsamples(D,':',':',':'));
	samples_of_interest = reshape(samples_of_interest,1,D.nsamples*D.ntrials);

	sm = D.ica.sm(D.ica.chan_inds,:); % For variance, only use the chan inds 
	tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';

	num_ics = D.ica.params.num_ics;

	%% COMPUTE METRICS 
	metrics = struct;
	
	% Mains frequency
	tc_fft = fft(demean(tc(:,samples_of_interest),2),[],2);
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
    kurt = kurtosis(tc(:,samples_of_interest),[],2);
    kurt = kurt(:);
    metrics.kurtosis.raw = kurt-3;
    metrics.kurtosis.log = log(kurt-3); 
    metrics.kurtosis.abs = abs(demean(boxcox1(kurt))); % make distribution more normal to better balance low/high kurtosis

	% Cardiac autocorrelation
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


	% Detect artefact channels related components
	if ~isempty(artefact_channels)

	    for artefact_chantype = unique(artefact_channels)
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

	% Add variance as default metric
	metrics.variance.value = nanvar(tc,[],2);
	metrics = orderfields(metrics);
	tc(:,~samples_of_interest) = NaN;


