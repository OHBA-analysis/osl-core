function badsections = osl_bad_sections(D,output_format)
% THIS FUNCTION WILL SOON BE REMOVED IN FAVOUR OF BADSAMPLES
% output format: 'logical' or 'indexed'
% AB 2013
%
% Changed by Giles Colclough 25 Nov 2013 to account for multiple trials. If
% there is more than one trial, D.events is a cell array of structures, one
% for each trial.

error(['Please replace use of this function with the SPM badsamples method... ', ...
      'e.g. badsamples(D,'':'','':'','':'')'])




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




    function badsections = get_bad_sections_from_events(Events, output_format)
        
        switch output_format
            case 'logical'
                if ~isempty(Events)
                    Events = Events(strcmp({Events.type},'BadEpoch'));
                    badsections = false(1,D.nsamples);
                    for ev = 1:numel(Events)
                        badsections = badsections | D.time >= Events(ev).time & D.time < (Events(ev).time+Events(ev).duration);
                    end
                else
                    badsections = false(size(D.time));
                end
                    
            case 'indexed'
                if ~isempty(Events)
                    Events = Events(strcmp({Events.type},'BadEpoch'));
                    
                    duration = [Events.duration];
                    time     = [Events.time];
                    if isempty(duration), duration = zeros(size(time)); end % catch case of duration being empty in each structure (it happens in some of the faces_subject_1 test data!)
                    
                    badsections = [time; time + duration]';
                else
                    badsections = [];
                end
        end%switch
        
    end%get_bad_sections_from_events

end
