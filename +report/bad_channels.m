function h = bad_channels(D)

	D = D.montage('switch',0);

	chaninds = D.badchannels;
	metric = std(D(:,:),[],2);
	modalities = unique(D.chantype);

	h = [];

	% MODALITY REPORTS
	for j = 1:length(modalities)

		h(end+1) = figure('name',sprintf('Bad Channels - modality ''%s'', %d',modalities{j}),'tag',sprintf('bad_channel_modality_%s',modalities{j}));
		this_modality = D.indchantype(modalities{j}); % All channel indices for the current modality
		this_modality_clean = setdiff(this_modality,chaninds);
		this_modality_bad = setdiff(this_modality,this_modality_clean);

		% Histogram
		subplot(4,1,1)
		[hs hsx]=hist(metric(this_modality),ceil(length(metric(this_modality))/10));
		bar(hsx,hs);
		title(sprintf('CHANNEL SUMMARY: %s', modalities{j}));
		a1=axis;
		xlabel('std')
		ylabel('counts')

		subplot(4,1,2);
		plot(metric(this_modality_clean),this_modality_clean,'og');
		if ~isempty(this_modality_bad)
			hold on
			plot(metric(this_modality_bad),this_modality_bad,'*r');
		end
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel('std')
		ylabel('channel index')

		subplot(4,1,3);
		[hs hsx]=hist(metric(this_modality_clean),ceil(length(metric(this_modality_clean))/10));
		bar(hsx,hs);
		a1=axis;
		title(sprintf('CHANNEL SUMMARY: %s (without bad channels)', modalities{j}));
		xlabel('std')
		ylabel('counts')

		subplot(4,1,4);
		plot(metric(this_modality_clean),this_modality_clean,'*g');
		a2=axis;
		axis([a1(1) a1(2) a2(3) a2(4) ]);
		xlabel('std')
		ylabel('channel index')
	
		badcolor     = [204     0     0] / 255;
		for j = 1:length(this_modality_bad)
			cl = D.chanlabels(this_modality_bad(j));
			cl = cl{1};
			unit = D.units(this_modality_bad(j));
			unit = unit{1};

			h(end+1) = figure('name',sprintf('Bad Channel - %s',cl),'tag',sprintf('bad_channel_%s',cl));
			pos = get(h(end),'Position');
			set(h(end),'Position',pos.*[1 1 1 0.5])
			plot(D.time,D(this_modality_bad(j),:),'Color',badcolor);
			title(sprintf('Bad Channel - %s',cl));
			xlabel('Time')
			ylabel(sprintf('Signal (%s)',unit));
		end

	end




