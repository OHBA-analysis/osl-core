function h = bad_trials(D,measure_fn)
	% Display summary for bad segments for each modality
	% Caution though - measure_fn or dummy_epoch_tsize could be different from what was used for automatic detection
	% However, recomputing the stats here means that this function is also able to coexist and work 
	% with manually identified bad segments

	assert(~strcmp(D.type,'continuous'),'Bad trial diagnostics only work for epoched MEEG objects');
	
	if nargin < 2 || isempty(measure_fn) 
		measure_fn = 'std';
	end
	
	D = D.montage('switch',0);
	candidate_modalities = {'EEG','MEG','MEGANY'};
	modalities = unique(D.chantype(D.indchantype(candidate_modalities)));

	h = [];
	for j = 1:length(modalities)
		chaninds = D.indchantype(modalities{j},'GOOD');
		data = D(chaninds,:,:);
		dat = feval(measure_fn,reshape(data,size(data,1)*size(data,2),size(data,3)),[],1); % Metric for each trial
		badtrials = ismember(1:D.ntrials,D.badtrials);

		h(end+1) = figure('name',sprintf('Bad trials - %s',modalities{j}),'tag',sprintf('bad_trials_modality_%s',modalities{j}));

		subplot(4,1,1)
		[hs hsx]=hist(dat,ceil(length(dat)/10));
		bar(hsx,hs);
		title(sprintf('TRIAL SUMMARY: %s', modalities{j}));
		a1=axis;
		xlabel(measure_fn)
		ylabel('counts')


		subplot(4,1,2);
		plot(dat(~badtrials),1:sum(~badtrials),'og');
		if any(badtrials)
			hold on
			plot(dat(badtrials),1:sum(badtrials),'*r');
		end
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel(measure_fn)
		ylabel('Epoch number')

		title(sprintf('Rejected %d of %d trials (%.2f%%)',sum(badtrials),length(badtrials),100*sum(badtrials)/length(badtrials)))

		subplot(4,1,3)
		[hs hsx]=hist(dat(~badtrials),ceil(length(dat(~badtrials))/10));
		bar(hsx,hs);
		a1=axis;
		xlabel(measure_fn)
		ylabel('counts')
        title(sprintf('TRIAL SUMMARY: %s (without bad channels)', modalities{j}));

		subplot(4,1,4);
		plot(dat(~badtrials),1:sum(~badtrials),'*g');
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel(measure_fn)
		ylabel('Epoch number')

	end



