function [D_epoched good_trial_starts] = osl_reject_bad_epoch_trials( S )

% [D_epoched good_trial_starts] = osl_reject_bad_epoch_trials( S )
%
% Adjust trials in D to reject those overlapping with bad epochs:
% Default is to assume S.events.type for bad epochs is 'BadEpoch'
% Note that S.events needs to contain the bad epochs.
%
% e.g.
%     S4=[];
%     S4.D_epoched=D_epoched;
%     S4.D_continuous=D_continuous;
%     S4.bad_event_type='BadEpoch';
%     epochinfo=osl_reject_bad_epoch_trials(S4);
%     figure;plot(D.time,good_trial_starts);
%
% MWW 2012

if ~isfield(S,'bad_event_type')
    S.bad_epoch_type='artefact_OPT';
end;

epochinfo=S.D_epoched.epochinfo;

eve = S.D_continuous.events; 

ind=0;start_bev=[];end_bev=[];
for ev = 1:numel(eve) % find bad epoch events
    if isfield(eve(ev),'type') && strcmp(eve(ev).type,S.bad_event_type) 
        ind=ind+1;
        start_bev(ind)=eve(ev).time;
        end_bev(ind)=eve(ev).time+eve(ev).duration;                               
    end
end

D_epoched=S.D_epoched;

if length(start_bev>0),% remove events overlapping with bad epochs:

    good_trial_starts=zeros(length(S.D_continuous.time),1);

    for ev = 1:size(epochinfo.trl,1) 
        if epochinfo.trl(ev,1) > 0 && epochinfo.trl(ev,2) <= length(S.D_continuous.time),

            st=S.D_continuous.time(epochinfo.trl(ev,1));
            en=S.D_continuous.time(epochinfo.trl(ev,2));
            keep=1;

            for bev=1:length(start_bev),

                if(st > start_bev(bev) && st < end_bev(bev))
                    keep=0;
                end;

                if(en > start_bev(bev) && en < end_bev(bev))
                    keep=0;
                end;

            end

            if ~keep,           
                D_epoched = badtrials(D_epoched, ev, 1); 
            else
                good_trial_starts(epochinfo.trl(ev,1))=1;            
            end;    
        else
            keep=0;
        end;
    end;      
else
    good_trial_starts=false(length(S.D_continuous.time),1);
    good_trial_starts(epochinfo.trl(:,1)) = true;
end;

end

