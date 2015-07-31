function glean_envelope(GLEAN)
% Runs the envelope stage of GLEAN

for session = 1:numel(GLEAN.data)
    
    if ~exist(GLEAN.data(session).enveloped,'file')
        
        % Make a temporary filename to copy raw data to
        [~,tempdata] = fileparts(tempname);
        tempdata = fullfile(fileparts(GLEAN.data(1).enveloped),[tempdata '.mat']);
        
        % Copy beamformed data to HMM directory with a temporary filename
        D = spm_eeg_load(GLEAN.data(session).beamformed);
        copy(D,tempdata);
        clear D
        
        % Set correct weights normalisation
        D = spm_eeg_load(tempdata);
        if montage(D,'getnumber') > 0
            if GLEAN.settings.envelope.weights_normalisation
                D = montage(D,'switch',2);
            else
                D = montage(D,'switch',1);
            end
            D.save;
        end
        
        % Parcellation & orthogonalisation
        if isstruct(GLEAN.settings.envelope.parcellation)
            S = [];
            S.D = tempdata;
            S.parcellation      = GLEAN.settings.envelope.parcellation.file;
            S.mask              = GLEAN.settings.envelope.parcellation.mask;
            S.orthogonalisation = GLEAN.settings.envelope.parcellation.orthogonalisation;
            S.method            = GLEAN.settings.envelope.parcellation.method;
            D = glean_parcellation(S);
            move(D,tempdata); % overwrite temporary file with enveloped data
        end
        
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
%        S.demean    = 1; % to ensure later normalisation/PCA is correct
        S.demean    = 0;
        S.prefix    = '';
        D = glean_hilbenv(S);
        
        % Rename file
        D = move(D,GLEAN.data(session).enveloped);
         
    end
    
end

end
