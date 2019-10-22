function hfig = trial_timings(D, epochinfo, plot_name_prefix)

% hfig = event_timings(D_epoched)
% or
% hfig = event_timings(D_continuous, epochinfo)
%
% Plots occurence of trials (whether they are marked as good or bad) and bad samples in unepoched, continuous
% time.
% Input D_epoched needs to be epoched data that contains an epochinfo field (normally this is set by osl_epoch)
%
% MWW

if nargin<3
    plot_name_prefix='';
end

if ischar(D)
    D=spm_eeg_load(D);
end

if nargin==1 || isempty(epochinfo)
    
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
    
    D_epoched=D;
    
    if D_epoched.fsample ~= D_continuous.fsample
        error('fsample must be same for D and D_continuous=D.epochinfo.D to plot trial_timings');
    end
else
    
    % 2 inputs
    if D.ntrials>1
        warning('Only works with two inputs, event_timings(D, epochinfo), if D is continuous data');
        return
    end
    D_continuous=D;
    
    % construct condlist from epochinfo
    condlist=unique(epochinfo.conditionlabels);
    
    D_epoched=[];
end
clear D;

% epochinfo.trl is [ begin_sample, end_sample, offset ] (in ms) of each 
% trial (e.g. computed by spm_eeg_definetrial)

% calculate vector of when each condition has a trial happening
condition_on=zeros(length(condlist),length(D_continuous.time));
condition_start=zeros(length(condlist),length(D_continuous.time));
condition_end=zeros(length(condlist),length(D_continuous.time));

yrange=[0, 0.5];

df=fliplr(linspace(yrange(1),yrange(2),length(condlist)+3));
for cc=1:length(condlist)            
            
    if isempty(D_epoched)
        for ee=1:size(epochinfo.trl,1)
            if epochinfo.trl(ee,1)>=1 && epochinfo.trl(ee,2)<=length(D_continuous.time)
                if strcmp(epochinfo.conditionlabels{ee},condlist{cc})
                    condition_on(cc,epochinfo.trl(ee,1):epochinfo.trl(ee,2))=df(cc+1);
                    condition_start(cc,epochinfo.trl(ee,1))=df(cc+1);
                    condition_end(cc,epochinfo.trl(ee,2))=df(cc+1);
                end
            end
        end
    else
        % note that the trials in epochinfo can be different to those in
        % D_epoched, due to trials that run off the start or the end of the
        % continuous data
        
        for ee=1:D_epoched.ntrials
            
            if strcmp(D_epoched.conditions{ee},condlist{cc})
                start=round(D_epoched.trialonset(ee)*D_epoched.fsample);

                condition_on(cc,start : start + length(D_epoched.time) - 1)=df(cc+1);
                condition_start(cc,start)=df(cc+1);
                condition_end(cc,start + length(D_epoched.time) - 1)=df(cc+1);
            end
            
        end
    end
end

% we now have all we need from D:
clear D_epoched;

% plot
hfig = figure('name',[plot_name_prefix 'Event timings'],'tag',[plot_name_prefix 'event_timings']);
hold on;
pp=[];

for cc=1:length(condlist)
    tmp=double(condition_on(cc,:));
    tmp(tmp==0)=nan;
    pp(cc)=plot(D_continuous.time, tmp, get_cols(cc), 'LineWidth',8);
    
    tmp=double(condition_start(cc,:));
    tmp(tmp==0)=nan;
    plot(D_continuous.time, tmp, [get_cols(cc) 'o'], 'MarkerSize',12, 'LineWidth',2);
    
    tmp=double(condition_end(cc,:));
    tmp(tmp==0)=nan;
    plot(D_continuous.time, tmp, [get_cols(cc) 'x'], 'MarkerSize',12, 'LineWidth',2);
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