function glean_envelope(GLEAN)
% Runs the envelope stage of GLEAN


for session = 1:numel(GLEAN.data)
    
    data     = GLEAN.data(session).beamformed;
    envelope = GLEAN.data(session).enveloped;
    
    if ~exist(envelope,'file')
        
        % Make a temporary filename to copy raw data to
        [~,tempdata] = fileparts(tempname);
        tempdata = fullfile(fileparts(envelope),[tempdata '.mat']);
        
        % Copy data to HMM directory with a temporary filename
        D = spm_eeg_load(data);
        copy(D,tempdata);
        clear D
        
        % Compute envelopes
        S               = [];
        S.D             = tempdata;
        S.fsample_new   = GLEAN.settings.envelope.fsample;
        S.logtrans      = GLEAN.settings.envelope.log;
        if isfield(GLEAN.settings.envelope,'freqbands')
            S.freqbands = GLEAN.settings.envelope.freqbands;
        else
            S.freqbands = [];
        end
        S.demean    = 0;
        S.prefix    = '';
        D = glean_hilbenv(S);
        
        % Rename file
        D = move(D,envelope);
        system(['rm ',strrep(tempdata,'.mat','.*at')])
        
    end
    
end

end
