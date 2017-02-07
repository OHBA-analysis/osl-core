function first_level_results = oat_stats_combine_grads( oat, first_level_results )

OSLDIR = getenv('OSLDIR');

first_level=oat.first_level;
% contrasts
contrast_list=oat.first_level.contrast;
for c=1:length(contrast_list),
    contrast_list{c}=contrast_list{c}(:);
end;

ntpts                       = length(first_level_results.times);
nfreqs                      = length(first_level_results.frequencies);

% get D and update path in case OAT dir has been moved
D = oat_get_sensordata( first_level_results );

% Check that grads haven't already been combined
if ~isempty(D.indchantype('MEGCMB'))
    warning('Grads have already been combined! skipping this step...');
    return
end

%% combine the planar gradiometers if in sensorspace
if strcmp(first_level_results.recon_method,'none')
    switch first_level.sensor_space_combine_planars
        case 'combine_cartesian'

            % get labels for planar 1, planar 2, and magnetometers
            lay_fn = fullfile(OSLDIR,'layouts','neuromag306planar3.lay');
            cfg=[];
            cfg.layout = lay_fn;
            layplan1 = ft_prepare_layout(cfg);

            lay_fn = fullfile(OSLDIR,'layouts','neuromag306planar2.lay');
            cfg=[];
            cfg.layout = lay_fn;
            layplan2 = ft_prepare_layout(cfg);

            lay_fn = fullfile(OSLDIR,'layouts','neuromag306mag.lay');
            cfg=[];
            cfg.layout = lay_fn;
            laygrad = ft_prepare_layout(cfg);

            lay_fn = fullfile(OSLDIR,'layouts','neuromag306cmb.lay');
            cfg=[];
            cfg.layout = lay_fn;
            laycmb = ft_prepare_layout(cfg);

            % get indices for planar 1, planar 2, and magnetometers
            [planar1labels,planar1inds_inD] = intersect(D.chanlabels,layplan1.label);
            [planar2labels,planar2inds_inD] = intersect(D.chanlabels,layplan2.label);
            [maglabels,maginds_inD]         = intersect(D.chanlabels,laygrad.label);

            magcope             = first_level_results.cope(maginds_inD,:,:,:);
            magstd              = first_level_results.stdcope(maginds_inD,:,:,:);
            if(isfield(first_level_results,'cope_by_state')),
                magstatecope        = first_level_results.cope_by_state(maginds_inD,:,:,:);
                magstatestd         = first_level_results.stdcope_by_state(maginds_inD,:,:,:);
            end;
            magpseudozstatvar   = first_level_results.pseudo_zstat_var(maginds_inD,:,:,:);

            % adopting the order of the planar1's, step through the
            % planar2's and find the corresponding pair.  If not
            % present, throw the planar1 data away.  If present,
            % combine the data from the two channels
            pairlist                = [];
            combinedCope            = nan(numel(planar1labels),size(first_level_results.cope,2),size(first_level_results.cope,3),size(first_level_results.cope,4));
            %combinedCopeByState     = nan(numel(planar1labels),size(first_level_results.cope,2),size(first_level_results.cope,3),size(first_level_results.cope,4));
            combinedStdCope         = nan(numel(planar1labels),size(first_level_results.cope,2),size(first_level_results.cope,3),size(first_level_results.cope,4));
            %combinedStdCopeByState  = nan(numel(planar1labels),size(first_level_results.cope,2),size(first_level_results.cope,3),size(first_level_results.cope,4));
            count = 0;
            clear newcmblabels
            for iFindPair = 1:numel(planar1inds_inD)
                correspondingInd = find(strncmpi(planar1labels{iFindPair},planar2labels,6));
                if ~isempty(correspondingInd)
                    count = count + 1;
                    pairlist(iFindPair,:) = [iFindPair,correspondingInd];

                    % combine the cope
                    planar1res_temp = first_level_results.cope(planar1inds_inD(iFindPair),:,:,:);
                    planar2res_temp = first_level_results.cope(planar2inds_inD(correspondingInd),:,:,:);
                    if (strcmp(oat.first_level.cope_type,'acope')) || (~(strcmp(oat.first_level.tf_method,'none'))) % if the copes are rectified
                        combinedCope(count,:,:,:) = (planar1res_temp+planar2res_temp)./2;
                        if count == 1
                            disp('Averaging the (A)COPEs for orthogonal planar gradiometers')
                        end
                    else
                        combinedCope(count,:,:,:) = sqrt( planar1res_temp.^2 + planar2res_temp.^2  );
                        if count == 1
                            warning('Gradiometers combined with ERF and COPE: gradiometer COPEs will be rectified');
                        end
                    end

                    % combine the stdcope
                    planar1res_temp = first_level_results.stdcope(planar1inds_inD(iFindPair),:,:,:);
                    planar2res_temp = first_level_results.stdcope(planar2inds_inD(correspondingInd),:,:,:);
                    combinedStdCope(count,:,:,:) = sqrt( planar1res_temp.^2 + planar2res_temp.^2  ); % same transformation here but for different reasons

                    % combine the cope_by_state
                    if(0),%if(NKglm>1)
                        planar1res_temp = first_level_results.cope_by_state(planar1inds_inD(iFindPair),:,:,:);
                        planar2res_temp = first_level_results.cope_by_state(planar2inds_inD(correspondingInd),:,:,:);

                        if (strcmp(oat.first_level.cope_type,'acope')) || (~(strcmp(oat.first_level.tf_method,'none'))) % if the copes are rectified
                            combinedCopeByState(count,:,:,:) = (planar1res_temp+planar2res_temp)./2;
                            if count == 1
                                disp('Averaging the (A)COPEs (by state) for orthogonal planar gradiometers')
                            end
                        else
                            combinedCopeByState(count,:,:,:) = sqrt( planar1res_temp.^2 + planar2res_temp.^2  );
                            %warning('Gradiometers combined with ERF and COPE: gradiometer COPEs will be rectified');
                        end

                        % combine the stdcope_by_state
                        planar1res_temp = first_level_results.stdcope_by_state(planar1inds_inD(iFindPair),:,:,:);
                        planar2res_temp = first_level_results.stdcope_by_state(planar2inds_inD(correspondingInd),:,:,:);
                    end;
                    %combinedStdCopeByState(count,:,:,:) = sqrt( planar1res_temp.^2 + planar2res_temp.^2  );

                    labind = find(strncmpi(planar1labels{iFindPair},laycmb.label,6));
                    newcmblabels{count} = laycmb.label{labind};
                else

                end
            end % for iFindPair = 1:numel(planar1inds_inD)

            combinedCope(isnan(combinedCope)) = [];
            combinedStdCope(isnan(combinedStdCope)) = [];
            %combinedCopeByState(isnan(combinedCopeByState)) = [];
            %combinedStdCopeByState(isnan(combinedStdCopeByState)) = [];

            % re-baseline
            if(sum(first_level.bc)>0),
                baseline_time_indices=first_level_results.times<first_level.baseline_timespan(2) & first_level_results.times>first_level.baseline_timespan(1);

                for c=1:length(contrast_list),
                    if(first_level.bc(c))
                        combinedCope(:,:,c,:) = combinedCope(:,:,c,:) - repmat( mean(combinedCope(:,baseline_time_indices,c,:),2) , [1 size(combinedCope,2) 1 1]);
                        %combinedCopeByState(:,:,c,:) = combinedCopeByState(:,:,c,:) - repmat( mean(combinedCopeByState(:,baseline_time_indices,c,:),2) , [1 size(combinedCopeByState,2) 1 1]);
                    end
                end
            end

            % write back into the first level cope
            first_level_results.cope                = [combinedCope; magcope];
            first_level_results.stdcope             = [combinedStdCope; magstd];
            if 0,%NKglm>1,
                first_level_results.cope_by_state       = [combinedCopeByState; magstatecope];
                first_level_results.stdcope_by_state    = [combinedStdCopeByState; magstatestd];
            end;
            first_level_results.chanind             = first_level_results.chanind(1:size(first_level_results.cope,1));
            first_level_results.mask_indices_in_source_recon = first_level_results.mask_indices_in_source_recon(1:size(first_level_results.cope,1));

            % modify the sensorspace D object stored in first_level_results
            if size(newcmblabels,1)~=102; newcmblabels=newcmblabels.';end
            if size(maglabels,1)~=102; maglabels=maglabels.';end

            % preserve Class channel label if necessary
            if ~isempty(strcmp(first_level_results.D_sensor_data.chanlabels,'Class'));
                classchanind=find(strcmp(first_level_results.D_sensor_data.chanlabels,'Class'));
            else
                classchanind = -1;
            end
            

            alllabels    = [newcmblabels; maglabels];
            first_level_results.D_sensor_data = chanlabels(first_level_results.D_sensor_data,1:numel(alllabels),alllabels);
            nExcessChans = numel(first_level_results.D_sensor_data.chanlabels) - numel(alllabels);
            first_level_results.D_sensor_data = chanlabels(first_level_results.D_sensor_data,numel(alllabels)+1:numel(alllabels)+nExcessChans,[]);
            megcmblab = {'MEGCMB'};
            megcmblab = repmat(megcmblab,[numel(newcmblabels),1]);
            first_level_results.D_sensor_data = chantype(first_level_results.D_sensor_data,1:numel(newcmblabels),megcmblab);
            megmaglab = {'MEGMAG'};
            magmaglab = repmat(megmaglab,[numel(maglabels),1]);
            first_level_results.D_sensor_data = chantype(first_level_results.D_sensor_data,numel(newcmblabels)+1:numel(alllabels),magmaglab);
            otherlab = {'Other'};
            otherlab = repmat(otherlab,[nExcessChans,1]);
            first_level_results.D_sensor_data = chantype(first_level_results.D_sensor_data,numel(alllabels)+1:numel(alllabels)+nExcessChans,otherlab);

            % Need to preserve this label if running sensorspace oat
            if classchanind > 0
                first_level_results.D_sensor_data = chanlabels(first_level_results.D_sensor_data,classchanind, 'Class');
            end

            first_level_results.D_sensor_data.save;

        case 'dont_combine'

    end % switch first_level.sensor_space_combine_planars

    chanindmeg = strmatch('MEG', first_level_results.D_sensor_data.chantype);
    first_level_results.chanind=chanindmeg;
    first_level_results.chanlabels=first_level_results.D_sensor_data.chanlabels(chanindmeg);
    first_level_results.chantype=first_level_results.D_sensor_data.chantype(chanindmeg);

end % if strcmp(source_recon_results.recon_method,'none')


end

