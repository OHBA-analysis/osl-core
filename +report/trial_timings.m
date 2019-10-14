function hfig = trial_timings(D, epochinfo)

% hfig = event_timings(D_epoched)
% or
% hfig = event_timings(D_continuous, epochinfo)
%
% Plots occurence of trials and bad samples in unepoched, continuous
% time.
% Input D_epoched needs to be epoched data that contains an epochinfo field (normally this is set by osl_epoch)
%
% MWW

if ischar(D)
    D=spm_eeg_load(D);
end

if nargin==1
    
    if D.ntrials==1
        warning('Only works with single input, event_timings(D), if D is epoched data');
        return
    end

    if ~isfield(D,'epochinfo')
        warning('Inputted epoched SPM MEEG object needs to have an epochinfo field (normally this is set by osl_epoch)');
        return
    end
    D_continuous=D.epochinfo.D;
    
    epochinfo=D.epochinfo;
    condlist=D.condlist;
    
else
    
    % 2 inputs
    if D.ntrials>1
        warning('Only works with two inputs, event_timings(D, epochinfo), if D is continuous data');
        return
    end
    D_continuous=D;
    
    % construct condlist from epochinfo
    condlist=unique(epochinfo.conditionlabels);
    
end

% epochinfo is the begin_sample, end_sample and offset (in ms) of each trial (e.g. computed by spm_eeg_definetrial)

% calculate vector of when each condition has a trial happening
condition_on=zeros(length(condlist),length(D_continuous.time));
condition_start=zeros(length(condlist),length(D_continuous.time));
condition_end=zeros(length(condlist),length(D_continuous.time));

yrange=[0, 0.5];

df=fliplr(linspace(yrange(1),yrange(2),length(condlist)+3));
for cc=1:length(condlist)            
            
    % note that the trials in epochinfo can be different to those in
    % D_epoched, due to trials that run of the start or the end of the
    % continuous data
    trial_onsets_trl=D.epochinfo.trl(:, 1)./D.fsample; % in secs

    for ee=1:D.ntrials
        % Find trial index in epochinfo by matching trial onsets
        ind=find(trial_onsets_trl==D.trialonset(ee));
        
        if strcmp(epochinfo.conditionlabels{ind},condlist{cc})
            condition_on(cc,epochinfo.trl(ind,1):epochinfo.trl(ind,2))=df(cc+1);
            condition_start(cc,epochinfo.trl(ind,1))=df(cc+1);
            condition_end(cc,epochinfo.trl(ind,2))=df(cc+1);
        end
    end

end

% we now have all we need from D:
clear D;

% plot
hfig = figure('name','Event timings','tag','event_timings');
hold on;
pp=[];

for cc=1:length(condlist)
    tmp=double(condition_on(cc,:));
    tmp(tmp==0)=nan;
    pp(cc)=plot(D_continuous.time, tmp, get_cols(cc), 'LineWidth',8);
    
    tmp=double(condition_start(cc,:));
    tmp(tmp==0)=nan;
    plot(D_continuous.time, tmp, ['ko'], 'MarkerSize',8, 'LineWidth',2);
    
    tmp=double(condition_end(cc,:));
    tmp(tmp==0)=nan;
    plot(D_continuous.time, tmp, ['kx'], 'MarkerSize',8, 'LineWidth',2);
end

tmp=double(~good_samples(D_continuous,D_continuous.indchantype('MEEG','GOOD')))*df(length(condlist)+2);
tmp(tmp==0)=nan;

pp(length(condlist)+1)=plot(D_continuous.time, tmp, 'k', 'LineWidth',8);

legend(pp,[condlist, {'Bad samples'}]);

set(gca,'YTick',[],'YTickLabel',[]);
title('Timings of trials and bad samples');
set(gca, 'YLim', yrange);
set(gca, 'XLim', [D_continuous.time(1) D_continuous.time(end)]);
plot4paper('Time (s)','');
set(hfig,'Position',[1 1 1500 400]);

end