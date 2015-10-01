function glean_results(GLEAN)
% Runs the results stage of GLEAN

model = load(GLEAN.model.model);

if isfield(GLEAN.envelope.settings,'freqbands')
    F = numel(GLEAN.envelope.settings.freqbands);
else
    F = 1;
end


for results_type = fieldnames(GLEAN.results.settings)'
    
    results = lower(char(results_type));
    
    switch results
        
        case 'pcorr'
            
            for subspace = fieldnames(GLEAN.results.pcorr)'
                
                switch char(subspace)
                    case 'voxel'
                        data = 'envelope';
                    case 'parcel'
                        data = 'subspace';
                end
                
                D = spm_eeg_load(GLEAN.(data).data{1});
                
                switch lower(char(fieldnames(GLEAN.model.settings)))
                    case 'hmm'
                        regressors = cell2mat(arrayfun(@(k) model.hmm.statepath==k,1:model.hmm.K,'UniformOutput',0));
                        session_maps = nan(D.nchannels,model.hmm.K,F,numel(GLEAN.data));
                    case 'ica'
                        regressors = model.ica.tICs';
                        session_maps = nan(D.nchannels,size(model.ica.tICs,1),F,numel(GLEAN.data));
                end
                
                
                for session = 1:numel(GLEAN.data)
                    
                    if F == 1
                        session_maps(:,:,1,session) = glean_regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),results);
                    else
                        session_maps(:,:,:,session) = glean_regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),results);
                    end
                    % Save the session specific maps
                    disp(['Saving partial correlation maps for session ' num2str(session)])
                    for f = 1:F
                        map = session_maps(:,:,f,session);
                        switch GLEAN.results.settings.(results).format
                            case 'mat'
                                save(GLEAN.results.(results).(char(subspace)).sessionmaps{session}{f},'map');
                            case 'nii'
                                save2nii(map,GLEAN.results.(results).(char(subspace)).sessionmaps{session}{f},char(subspace))
                        end
                    end
                    
                end
                group_maps = nanmean(session_maps,4);
                
                
                % Save the group averaged maps
                disp('Saving group partial correlation map')
                for f = 1:F
                    map = group_maps(:,:,f);
                    switch GLEAN.results.settings.(results).format
                        case 'mat'
                            save(GLEAN.results.(results).(char(subspace)).groupmaps{f},'map');
                        case 'nii'
                            save2nii(map,GLEAN.results.(results).(char(subspace)).groupmaps{f},char(subspace))
                    end
                end
            end
            
        case 'connectivity_profile'
            if F > 1
                error('Not yet supported for multiband HMM');
            end
            
            group_maps = glean_connectivityprofile(model.hmm);
            
            % Save the group averaged maps
            disp('Saving group connectivity_profile map')
            for f = 1:F
                map = group_maps(:,:,f);
                switch GLEAN.results.settings.(results).format
                    case 'mat'
                        save(GLEAN.results.(results).groupmaps{f},'map');
                    case 'nii'
                        save2nii(map,GLEAN.results.(results).groupmaps{f},'parcel')
                end
            end
            
            
    end
    
end



    function save2nii(map,fname,space)
        % Have to work out what spatial basis set we're in:
        % pre-envelope parcellation (orthogonalisation) - use parcellation as mask
        % post-envelope parcellation or pca - use full voxelwise mask for pcorr,
        %                                   - use parcellation for connectivity profile
        
        switch space
            case 'voxel'
                writenii(map,fname,GLEAN.results.settings.(results).mask);
            case 'parcel'
                map = parcellation2map(map,GLEAN.subspace.settings.parcellation.file,GLEAN.subspace.settings.parcellation.mask);
                writenii(map,fname,GLEAN.results.settings.(results).mask);
        end

    end

end