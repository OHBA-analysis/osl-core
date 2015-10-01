function glean_envelope(GLEAN)
% Runs the envelope stage of GLEAN

for session = 1:numel(GLEAN.data)
    
    if ~exist(GLEAN.envelope.data{session},'file')
        
        % Make a temporary filename to copy raw data to
        [~,tempdata] = fileparts(tempname);
        tempdata = fullfile(fileparts(GLEAN.envelope.data{session}),[tempdata '.mat']);
        
        % Copy data to HMM directory with a temporary filename
        D = spm_eeg_load(GLEAN.data{session});
        copy(D,tempdata);
        clear D
        
        % Compute envelopes
        S               = [];
        S.D             = tempdata;
        S.fsample_new   = GLEAN.envelope.settings.fsample;
        S.logtrans      = GLEAN.envelope.settings.log;
        if isfield(GLEAN.envelope.settings,'freqbands')
            S.freqbands = GLEAN.envelope.settings.freqbands;
        else
            S.freqbands = [];
        end
        S.demean    = 0;
        S.prefix    = '';
        D = glean_hilbenv(S);
        
        % Rename file
        move(D,GLEAN.envelope.data{session});
        
        % Tidy up
        system(['rm ',strrep(tempdata,'.mat','.*at')])

    end
        
end

end
