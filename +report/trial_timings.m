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

else
    % 2 inputs
    if D.ntrials>1
        warning('Only works with two inputs, event_timings(D, epochinfo), if D is continuous data');
        return
    end
    D_continuous=D;
end


% epochinfo is the begin_sample, end_sample and offset (in ms) of each trial (e.g. computed by spm_eeg_definetrial)

% calculate vectro of when each condition has a trial happening
condition_on=zeros(length(D.condlist),length(D_continuous.time));

yrange=[0, 0.5];

df=linspace(yrange(1),yrange(2),length(D.condlist)+3);
for cc=1:length(D.condlist)
    
    for ee=1:length(epochinfo.conditionlabels)
        if strcmp(epochinfo.conditionlabels{ee},D.condlist{cc})
            condition_on(cc,epochinfo.trl(ee,1):epochinfo.trl(ee,2))=df(cc+1);
        end
    end

end

% plot
hfig = figure('name','Event timings','tag','event_timings');
hold on;

for cc=1:length(D.condlist)
    tmp=double(condition_on(cc,:));
    tmp(tmp==0)=nan;
    plot(D_continuous.time, tmp, get_cols(cc), 'LineWidth',8);
end

tmp=double(~good_samples(D_continuous,D_continuous.indchantype('MEEG','GOOD')))*df(length(D.condlist)+2);
tmp(tmp==0)=nan;

plot(D_continuous.time, tmp, 'k', 'LineWidth',8);

legend([D.condlist, {'Bad samples'}]);

set(gca,'YTick',[],'YTickLabel',[]);
title('Timings of trials and bad samples');
set(gca, 'YLim', yrange);
plot4paper('Time (s)','');

end