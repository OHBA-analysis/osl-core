function h = bad_segments(D)
	% Display summary for bad segments for each modality

	D = D.montage('switch',0);
	measure_fn = 'std';
	dummy_epoch_tsize = 1;
	modalities = unique(D.chantype);
	dummy_trialsize = dummy_epoch_tsize*D.fsample;
	dummy_ntrials = floor(D.nsamples/(dummy_epoch_tsize*D.fsample));
	convert_to_trial = @(x) reshape(x(:,1:dummy_trialsize*dummy_ntrials),size(x,1),dummy_trialsize,dummy_ntrials); % This function 'epochs' a matrix from chans x continuous_time to chans x trial_time x trials

	h = [];
	for j = 1:length(modalities)
		chaninds = D.indchantype(modalities{j},'GOOD');
		data = convert_to_trial(D(chaninds,:));
		badtrials = squeeze(any(convert_to_trial(~good_samples(D,chaninds)),2));
		dat = feval(measure_fn,reshape(data,size(data,1)*size(data,2),size(data,3)),[],1); % Metric for each virtual trial

		h(end+1) = figure('name',sprintf('Bad Segments - modality ''%s'', %d',modalities{j}),'tag',sprintf('bad_segments_modality_%s',modalities{j}));

		subplot(4,1,1)
		[hs hsx]=hist(dat,ceil(length(dat)/10));
		bar(hsx,hs);
		title(sprintf('VIRTUAL EPOCH (%.2fs) SUMMARY: %s', dummy_epoch_tsize, modalities{j}));
		a1=axis;
		xlabel('std')
		ylabel('counts')


		subplot(4,1,2);
		plot(dat(~badtrials),1:sum(~badtrials),'og');
		if any(badtrials)
			hold on
			plot(dat(badtrials),1:sum(badtrials),'*r');
		end
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel('std')
		ylabel('Epoch number')

		title(sprintf('Rejected %d of %d trials (%.2f%%)',sum(badtrials),length(badtrials),100*sum(good_samples(D,chaninds))/D.nsamples))

		subplot(4,1,3)
		[hs hsx]=hist(dat(~badtrials),ceil(length(dat(~badtrials))/10));
		bar(hsx,hs);
		a1=axis;
		xlabel('std')
		ylabel('counts')


		subplot(4,1,4);
		plot(dat(~badtrials),1:sum(~badtrials),'*g');
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel('std')
		ylabel('Epoch number')

	end



