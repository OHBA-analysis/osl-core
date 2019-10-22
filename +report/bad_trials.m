function h = bad_trials(D,measure_fn,modalities,plot_name_prefix)
    
% h = bad_trials(D)
% h = bad_trials(D,measure_fn)
% h = bad_trials(D,measure_fn,modalities)
% h = bad_trials(D,measure_fn,modalities,plot_name_prefix)
%
% Display summary of epoched SPM MEG object D for bad trials for each modality

    assert(~strcmp(D.type,'continuous'),'Bad trial diagnostics only work for epoched MEEG objects');
    
    if nargin < 2 || isempty(measure_fn) 
        measure_fn = 'std';
    end
    
    if nargin<3
        candidate_modalities = {'EEG','MEG','MEGANY'};
        modalities = unique(D.chantype(D.indchantype(candidate_modalities)));
    end
    
    if nargin<4
        plot_name_prefix='';
    end
    
    D = D.montage('switch',0);
    
    h = [];
    for j = 1:length(modalities)
        chaninds = D.indchantype(modalities{j},'GOOD');
        data = D(chaninds,:,:);
        dat = feval(measure_fn,reshape(data,size(data,1)*size(data,2),size(data,3)),[],1); % Metric for each trial
        badtrials = ismember(1:D.ntrials,D.badtrials);

        h(end+1) = figure('name',sprintf([plot_name_prefix 'Bad trials - ' modalities{j}]),'tag',sprintf([plot_name_prefix 'bad_trials_modality_' modalities{j}]));

        subplot(4,1,1)
        [hs hsx]=hist(dat,ceil(length(dat)/10));
        bar(hsx,hs);
        title(sprintf('TRIAL SUMMARY: %s', modalities{j}));
        a1=axis;

        subplot(4,1,2);
        plot(dat(~badtrials),find(~badtrials),'og');
        if any(badtrials)
            hold on
            plot(dat(badtrials),find(badtrials),'*r');
        end
        a2=axis;
        axis([a1(1) a1(2) a2(3) a2(4) ]);
        ylabel('Trial index');
    
        title(sprintf('Rejected %d of %d trials (%.2f%%)',sum(badtrials),length(badtrials),100*sum(badtrials)/length(badtrials)))

        subplot(4,1,3)
        [hs hsx]=hist(dat(~badtrials),ceil(length(dat(~badtrials))/10));
        bar(hsx,hs);
        a1=axis;
        
        title(sprintf('TRIAL SUMMARY: %s (without bad trials)', modalities{j}));

        subplot(4,1,4);
        plot(dat(~badtrials),find(~badtrials),'*g');
        a2=axis;
        axis([a1(1) a1(2) a2(3) a2(4) ]);
        xlabel(measure_fn);
        ylabel('Trial index');
    
        set(h(end),'Position',[1 1 500 700]);

    end



