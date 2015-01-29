function [ D_epoched goodtrials ] = osl_epoch( S )

% [ D_epoched goodtrials ] = osl_epoch( S )
%
% epochs data using S.epochinfo and 
% (if reject_bad_epochs flag set) rejects any epochs overlapping with bad
% events, where bad events are indicated by S.D.event.type being the same
% as S.bad_epoch_type.
%
% Returns epoched data, D_epoched. The S.epochinfo is stored inside
% D_epoched.
%
% e.g.
% S=[];
% S.D=D; % continuous data
% S.epochinfo=epochinfo;
% [ D_epoched goodtrials ] = osl_epoch( S )
%
% MWW 2012

if ~isfield(S,'bad_event_type')
    S.bad_event_type='artefact_OSL';
end;

if ~isfield(S,'reject_bad_epochs'),
    S.reject_bad_epochs=1;
end;

S3=[];
S3.D = S.D;    
S3.save=0;
S3.reviewtrials=0;
S3.bc=0;     
S3.trl=S.epochinfo.trl;
S3.conditionlabels=S.epochinfo.conditionlabels;

D_epoched = spm_eeg_epochs(S3);

D_epoched.epochinfo=S.epochinfo; % store epoch info inside D object  
D_epoched.epochinfo.time_continuous=S.D.time;

if(S.reject_bad_epochs)
    S4=[];
    S4.D_epoched=D_epoched;
    S4.D_continuous=S.D;
    S4.bad_event_type=S.bad_event_type;
    [D_epoched goodtrials]=osl_reject_bad_epoch_trials(S4);
else 
    goodtrials=[];
end;

D_epoched.save;

end

