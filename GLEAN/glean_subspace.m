function glean_subspace(GLEAN)
% Runs the subspace stage of GLEAN

if ~all(cellfun(@exist,{GLEAN.data.subspace}))
      
    method = lower(char(intersect(fieldnames(GLEAN.settings.subspace),{'pca','parcellation','voxel'})));
        
    % This should probs be a temporary file
    for session = 1:numel(GLEAN.data)
        % Copy envelope data to subspace directory
        D = spm_eeg_load(GLEAN.data(session).enveloped);
        D = copy(D,prefix(GLEAN.data(session).subspace,'tmp'));
        
        % Apply normalisation
        S               = [];
        S.D             = prefix(GLEAN.data(session).subspace,'tmp');
        S.normalisation = GLEAN.settings.subspace.normalisation;
        glean_normalise(S);

        %D(:,:,:,:) = D(:,:,:,:) ./ repmat(D.normalisation,1,D.nfrequencies,D.nsamples,D.ntrials);
        %D.save;

    end
    
    
    switch method
        
        case 'voxel'
            
            % Just return the identity matrix:
            D = spm_eeg_load(GLEAN.data(1).subspace);
            M = sparse(eye(D.nchannels));
            
        case 'pca'
            
            C = osl_groupcov(prefix({GLEAN.data.subspace},'tmp'));
            pcadim = min(GLEAN.settings.subspace.pca.dimensionality,D.nchannels);
            [allsvd,M] = eigdec(C,pcadim);
            
            if GLEAN.settings.subspace.pca.whiten
                M = diag(1./sqrt(allsvd)) * M';
            else
                M = M';
            end
            
        case 'parcellation'
            
            % Parcellation - pass in as a P*V matrix or 1*V vector
            parcellationFile = GLEAN.settings.subspace.parcellation.file;
            parcellationMask = GLEAN.settings.subspace.parcellation.mask;
            if ~isempty(strfind(parcellationFile,'.nii'))
                parcellation = readnii(parcellationFile,parcellationMask);
            elseif ~isempty(strfind(parcellationFile,'.mat'))
                parcellation = load(parcellationFile);
                if length(fieldnames(parcellation)) == 1
                    parcellation = parcellation.(char(fieldnames(parcellation)));
                else
                    error('.mat file should contain only one variable')
                end
            end
            
            switch GLEAN.settings.subspace.parcellation.method
                case 'spatialbasis'
                    M = parcellation';
                    M(isnan(M)) = 0;
                otherwise
                    error('method not supported');
            end
            
            
        otherwise
            error('I don''t know what to do!')
            
    end
    
    
    % Apply spatial basis and write output files
    for session = 1:numel(GLEAN.data)
        
        D = spm_eeg_load(prefix(GLEAN.data(session).subspace,'tmp'));
        
        montnew             = [];
        montnew.name        = 'spatialbasis';
        montnew.labelnew    = arrayfun(@(x) strcat(method,num2str(x)),1:size(M,1),'uniformoutput',0)';
        montnew.labelorg    = D.chanlabels;
        montnew.tra         = M;     
        
        S2 = [];
        S2.D            = prefix(GLEAN.data(session).subspace,'tmp');
        S2.montage      = montnew;
        S2.keepsensors  = false;
        S2.keepothers   = false;
        S2.mode         = 'write';
        
        D = spm_eeg_montage(S2);
        D.save;
        
        move(D,GLEAN.data(session).subspace)
        unix(['rm ' strrep(prefix(GLEAN.data(session).subspace,'tmp'),'.mat','.*at')]);
    end
    
end





        

