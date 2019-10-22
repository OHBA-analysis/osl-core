function h = bad_segments(D,S)

% h = bad_segments(D,S)
%
% S.dummy_epoch_tsize
% S.modalities
% S.measure_fn
%
% Display summary of continuous SPM MEG object D for bad segments for each modality
% Caution though - measure_fn or dummy_epoch_tsize could be different from what was used for automatic detection
% However, recomputing the stats here means that this function is also able to coexist and work 
% with manually identified bad segments

assert(strcmp(D.type,'continuous'),'Bad segment diagnostics only work for continuous MEEG objects');

if nargin < 2
    S=[];
end

try plot_name_prefix=S.plot_name_prefix; catch, plot_name_prefix = ''; end
try dummy_epoch_tsize=S.dummy_epoch_tsize; catch, dummy_epoch_tsize = 1; end
try modalities=S.modalities; 
catch    
    candidate_modalities = {'EEG','MEG','MEGANY'};
    modalities = unique(D.chantype(D.indchantype(candidate_modalities)));
end
try measure_fn=S.measure_fn; catch, measure_fn = 'std'; end
try do_virtual_epoch_plots=S.do_virtual_epoch_plots; catch, do_virtual_epoch_plots = false; end

D = D.montage('switch',0);

dummy_trialsize = dummy_epoch_tsize*D.fsample;
dummy_ntrials = floor(D.nsamples/(dummy_epoch_tsize*D.fsample));
convert_to_trial = @(x) reshape(x(:,1:dummy_trialsize*dummy_ntrials),size(x,1),dummy_trialsize,dummy_ntrials); % This function 'epochs' a matrix from chans x continuous_time to chans x trial_time x trials

h = [];
for j = 1:length(modalities)
    if do_virtual_epoch_plots
        chaninds = D.indchantype(modalities{j},'GOOD');
        data = convert_to_trial(D(chaninds,:));
        badtrials = squeeze(any(convert_to_trial(~good_samples(D,chaninds)),2));
        dat = feval(measure_fn,reshape(data,size(data,1)*size(data,2),size(data,3)),[],1); % Metric for each virtual trial

    
        h(end+1) = figure('name',sprintf([plot_name_prefix 'Bad Segments - ' modalities{j}]),'tag',sprintf([plot_name_prefix 'bad_segments_modality_' modalities{j}]));

        subplot(4,1,1)
        [hs hsx]=hist(dat,ceil(length(dat)/10));
        bar(hsx,hs);
        title(sprintf('VIRTUAL EPOCH (%.2fs) SUMMARY: %s', dummy_epoch_tsize, modalities{j}));
        a1=axis;

        subplot(4,1,2);
        plot(dat(~badtrials),find(~badtrials),'og');
        if any(badtrials)
            hold on
            plot(dat(badtrials),find(badtrials),'*r');
        end
        a2=axis;
        axis([a1(1) a1(2) a2(3) a2(4) ]);
        ylabel('Epoch number');

        title(sprintf('Rejected %d of %d trials (%.2f%%)',sum(badtrials),length(badtrials),100*sum(good_samples(D,chaninds))/D.nsamples))

        subplot(4,1,3)
        [hs hsx]=hist(dat(~badtrials),ceil(length(dat(~badtrials))/10));
        bar(hsx,hs);
        a1=axis;
        title(sprintf('VIRTUAL EPOCH (%.2fs) SUMMARY: %s (without bad segments)', dummy_epoch_tsize, modalities{j}));

        subplot(4,1,4);
        plot(dat(~badtrials),find(~badtrials),'*g');
        a2=axis;
        axis([a1(1) a1(2) a2(3) a2(4) ]);
        xlabel(measure_fn);
        ylabel('Epoch number');

        set(h(end),'Position',[1 1 500 700]);

    end
    
    h(end+1) = figure('name',sprintf([plot_name_prefix 'Bad segment timings - ' modalities{j}]),'tag',sprintf([plot_name_prefix 'bad_segments_timings_modality_' modalities{j}]));

    tmp=double(~good_samples(D,D.indchantype(modalities{j},'GOOD'))*0.05);
    tmp(tmp==0)=nan;

    plot(D.time, tmp, 'k', 'LineWidth', 16);

    set(gca,'YTick',[],'YTickLabel',[]);
    set(gca, 'YLim', [0 0.1]);
    set(gca, 'XLim', [D.time(1) D.time(end)]);
    set(h(end),'Position',[1 1 1500 150]);
    plot4paper('Time (s)','');
    
end



