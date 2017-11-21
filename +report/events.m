function hfig = events(D)
    % Display summary information about events
    
    hfig = figure('name','Events','tag','events_hist');
    
    ev = events(D,1);
    if ~isempty(ev)
        ev = ev(cellfun(@(x) isempty(strmatch('artefact',x)),{ev.type})); % Select all non-artefact events
    end
    
    if isempty(ev)
        title('No events found');
        return
    end
    
    vals = [ev.value];
    counts = accumarray(vals(:),1)/2; % Divide by 2 because both up and down trigger
    vals = find(counts);

    bar(1:length(vals),counts(vals));
    hold on
    for j = 1:length(vals)
        text(j,counts(vals(j)),num2str(counts(vals(j))),'HorizontalAlignment','center','VerticalAlignment','bottom','Color','r')
    end

    set(gca,'XTick',1:length(vals),'XTickLabel',vals);
    fig_names{1}='events_hist';
    title('Non-artefact event counts corrected for button presses');
    xlabel('Trigger codes');
    ylabel('# of triggers');
    set(gca,'YLim',[0 ceil(max(counts)*1.1)]); % Increase Y limit by 10%

end