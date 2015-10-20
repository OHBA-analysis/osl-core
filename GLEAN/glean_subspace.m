function glean_subspace(GLEAN)
% Runs the subspace stage of GLEAN.
%
% GLEAN = GLEAN_SUBSPACE(GLEAN)
%
% Adam Baker 2015

pretty_string('RUNNING SUBSPACE STAGE')

method = lower(char(intersect(fieldnames(GLEAN.subspace.settings),{'pca','parcellation','voxel'})));

% Check if envelope file exists and whether or not to overwrite
files_exists = all(cellfun(@(file) exist(file,'file')==2, GLEAN.subspace.data));
overwrite   = GLEAN.subspace.settings.overwrite == 1;
if files_exists
    if overwrite
        msg = ['Overwriting existing subspace files in: \n' fullfile(GLEAN.subspace.dir,'data') '\n'];
        run_stage = true;
    else
        msg = ['Using existing subspace files in: \n' fullfile(GLEAN.subspace.dir,'data') '\n'];
        run_stage = false;
    end
else
    msg = ['Creating new subspace files in: \n' fullfile(GLEAN.subspace.dir,'data') '\n'];
    run_stage = true;
end
fprintf(msg);


if run_stage
      
    % --- Copy data to subspace directory ---
    for session = 1:numel(GLEAN.data)
        switch method
            case {'voxel','pca'}
                D = spm_eeg_load(GLEAN.envelope.data{session});
            case 'parcellation'
                D = spm_eeg_load(GLEAN.data{session});
        end
        
        D = copy(D,GLEAN.subspace.data{session});
    end
    
    % --- Apply normalisation ---
    for session = 1:numel(GLEAN.data)
        switch GLEAN.subspace.settings.normalisation
            
            case {'voxel','global'}
                D = spm_eeg_load(GLEAN.subspace.data{session});
                stdev = sqrt(glean_variance(D));
                if strcmp(GLEAN.subspace.settings.normalisation,'global')
                    stdev = mean(stdev);
                end
                M = montage(D,'getmontage');
                M.tra = diag(1./stdev)*M.tra;  
                % Remove unused montages:
                D = montage(D,'remove',1:montage(D,'getnumber'));
                D = montage(D,'add',M);
                D.save;
                
            case 'none'
                % Do nothing
                
            otherwise
                error('Invalid normalisation')
                
        end 
    end
    
    % --- Compute subspace ---
    switch method
        
        case 'voxel'
            % Do nothing
            
        case 'pca'
            C = glean_groupcov(GLEAN.subspace.data);
            pcadim = min(GLEAN.subspace.settings.pca.dimensionality,D.nchannels);
            [allsvd,M] = eigdec(C,pcadim);
            if GLEAN.subspace.settings.pca.whiten
                M = diag(1./sqrt(allsvd)) * M';
            else
                M = M';
            end        
            % Apply spatial basis and write output files
            for session = 1:numel(GLEAN.data)
                D = spm_eeg_load(GLEAN.subspace.data{session});
                S                  = [];
                S.D                = GLEAN.subspace.data{session};
                S.montage          = [];
                S.montage.name     = method;
                S.montage.labelnew = arrayfun(@(x) strcat(method,num2str(x)),1:size(M,1),'uniformoutput',0)';
                S.montage.labelorg = D.chanlabels;
                S.montage.tra      = M;
                S.keepsensors      = false;
                S.keepothers       = false;
                S.mode             = 'write';
                D = spm_eeg_montage(S);
                move(D,GLEAN.subspace.data{session});
            end   
            
        case 'parcellation'
            for session = 1:numel(GLEAN.data)   
                % Compute parcellation:
                S                   = [];
                S.D                 = GLEAN.subspace.data{session};
                S.parcellation      = GLEAN.subspace.settings.parcellation.file;
                S.mask              = GLEAN.envelope.settings.mask;
                S.orthogonalisation = GLEAN.subspace.settings.parcellation.orthogonalisation;
                S.method            = GLEAN.subspace.settings.parcellation.method;
                glean_parcellation(S);
                % Compute envelopes:
                S               = [];
                S.D             = GLEAN.subspace.data{session};
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
                move(D,GLEAN.subspace.data{session});
            end
            
        otherwise
            error('I don''t know what to do!')
            
    end
    
end
