function plot_bad_components(D,tc,varargin)


	arg.addParameter('do_mains',false); 
	arg.addParameter('mains_frequency',50); 
	arg.addParameter('do_kurt',false); 
	arg.addParameter('do_cardiac',false); 



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
	if IDENT_PARAMS123.do_mains
	    spec_n = normalise(spec,2);
	    inds = freq_ax > IDENT_PARAMS123.mains_frequency - 1 & freq_ax < IDENT_PARAMS123.mains_frequency + 1;
	    metrics.mains.value = max(spec_n(:,inds),[],2);
	end

	% Kurtosis 
	if IDENT_PARAMS123.do_kurt
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
	if ~isempty(S.artefact_channels)

	    for artefact_chantype = unique(S.artefact_channels)
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

