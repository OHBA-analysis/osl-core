function badsections = AB_get_bad_sections(D,output_format)
% output format: 'logical' or 'indexed'

% AB 2013

% Changed by Giles Colclough 25 Nov 2013 to account for multiple trials. If
% there is more than one trial, D.events is a cell array of structures, one
% for each trial. 

if nargin<2
  output_format = 'indexed';
end

if D.ntrials > 1,
    badsections = [];
    Events = events(D, ':');
    for iTrial = 1:D.ntrials,
        badsections = [badsections; ...
                       get_bad_sections_from_events(Events{iTrial}, ...
                                                    output_format)];       %#ok<AGROW>
    end%if
    
elseif D.ntrials == 1,
    badsections = get_bad_sections_from_events(events(D,1), output_format);
    
else
    error('What''s up with the number of trials?\n');
end%if

end%AB_get_bad_sections



function badsections = get_bad_sections_from_events(Events, output_format)
if ~isempty(Events)
  Events = Events(strcmp({Events.type},'BadEpoch'));
  
  switch output_format
    case 'logical'
      badsections = false(1,D.nsamples);
      for ev = 1:numel(Events)
        badsections = badsections | D.time >= Events(ev).time & D.time < (Events(ev).time+Events(ev).duration);
      end
    case 'indexed'
        duration = [Events.duration]; 
        time     = [Events.time];
        if isempty(duration), duration = zeros(size(time)); end % catch case of duration being empty in each structure (it happens in some of the faces_subject_1 test data!)
        
        badsections = [time; time + duration]';
      
  end%switch

else
  badsections = [];
  
end%if

end%get_bad_sections_from_events