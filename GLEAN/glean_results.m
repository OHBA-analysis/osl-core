function glean_results(GLEAN)
% Runs the results stage of GLEAN.
%
% GLEAN_RESULTS(GLEAN)
%
% Adam Baker 2015


pretty_string('RUNNING RESULTS STAGE')

model = load(GLEAN.model.model);

if isfield(GLEAN.envelope.settings,'freqbands')
    F = numel(GLEAN.envelope.settings.freqbands);
else
    F = 1;
end

for results_type = setdiff(fieldnames(GLEAN.results.settings),'dir')'
    
    results = lower(char(results_type));
    
    switch results
        
        case 'pcorr' % Partial correlation with inferred time courses
            
            for subspace = fieldnames(GLEAN.results.pcorr)'
                
                switch char(subspace)
                    case 'voxel'
                        data = 'envelope';
                    case 'parcel'
                        data = 'subspace';
                end
                
                D = spm_eeg_load(GLEAN.(data).data{1});
                
                % Regressors are the state time courses (HMM) or independent components (ICA)
                switch char(intersect(lower(fieldnames(GLEAN.model.settings)),{'hmm','ica'}));
                    case 'hmm'
                        regressors = cell2mat(arrayfun(@(k) model.hmm.statepath==k,1:model.hmm.K,'UniformOutput',0));
                        session_maps = nan(D.nchannels,model.hmm.K,F,numel(GLEAN.data));
                    case 'ica'
                        regressors = model.ica.tICs';
                        session_maps = nan(D.nchannels,size(model.ica.tICs,1),F,numel(GLEAN.data));
                end
                
                % Compute partial correlation map for each subject
                for session = 1:numel(GLEAN.data)
                    
                    if F == 1
                        session_maps(:,:,1,session) = glean_regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),results);
                    else
                        session_maps(:,:,:,session) = glean_regress(GLEAN.(data).data{session},regressors(model.subIndx==session,:),results);
                    end
                    % Save the session specific maps
                    disp(['Saving ' char(subspace) ' partial correlation maps for session ' num2str(session)])

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
                % Compute group map as the average of the session maps
                group_maps = nanmean(session_maps,4);
                
                % Save the group averaged maps
                disp(['Saving ' char(subspace) ' group partial correlation map'])
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
                        save(GLEAN.results.(results),'map');
                    case 'nii'
                        save2nii(map,GLEAN.results.(results),'parcel')
                end
            end
            
    end
    
end



    function save2nii(map,fname,space)
    % Save a matrix MAP to a nifti file with filename FNAME using a mask
    % appropriate for the SPACE the map is in (voxelwise or parcelwise)
        switch space
            case 'voxel'
                writenii(map,fname,GLEAN.envelope.settings.mask);
            case 'parcel'
                map = parcellation2map(map,GLEAN.subspace.settings.parcellation.file,GLEAN.envelope.settings.mask);
                writenii(map,fname,GLEAN.envelope.settings.mask);
        end

    end

end