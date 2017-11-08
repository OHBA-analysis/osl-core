function hfig = events(D)
    % Display summary information about events
    
    ev = events(D,1);
    ev = ev(cellfun(@(x) isempty(strmatch('artefact',x)),{ev.type})); % Select all non-artefact events
    vals = [ev.value];
    bins=unique(vals); % HL mod 1.1
    bins=[bins(1)-1 bins bins(end)+1]; % HL mod 1.1
    h=hist(vals,bins); % HL mod 1.1
    
    hfig = figure('name','Events','tag','events_hist');
    bar(bins,h/2); % HL mod 1.1, correct by 2 due to up and down triggers
    fig_names{1}='events_hist';
    title('Histogram of non-artefact events corrected for button presses');
    xlabel('Trigger codes');
    ylabel('# of triggers');
    set(gca,'YLim',[0 ceil(max(h/2)*1.1)]); % Increase Y limit by 10%

end