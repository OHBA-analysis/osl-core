function [ results_fnames ] = osl_run_subject_level( oat )

% [ results_fnames ] = osl_run_subject_level( oat )
%
% takes in an OAT, which needs to be setup by calling oat=osl_check_oat(S), struct
% and runs subject level session-wise fixed effects averaging
%
% This function should normally be called using osl_run_oat(oat);
%
% MWW 2011

OSLDIR = getenv('OSLDIR');

subject_level_results=[];

just_mean=1;

if max(oat.subject_level.subjects_to_do)>length(oat.subject_level.session_index_list),
    error('oat.subject_level.subjects_to_do and oat.subject_level.session_index_list have incompatible numbers of subjects');
end;

% do each subject separately
for subi=1:length(oat.subject_level.subjects_to_do),

    sub=oat.subject_level.subjects_to_do(subi);

    disp(['Doing fixed effects averaging over sessions for subject ' num2str(sub)]);

    num_sessions=length(oat.subject_level.session_index_list{sub});

    for sesi=1:num_sessions,

        ses=oat.subject_level.session_index_list{sub}(sesi);

        % load in single session stats
        first_level_results=oat_load_results(oat, oat.first_level.results_fnames{ses});

        if(sesi==1)
            sc=[size(first_level_results.cope) 1 1];
            cope=zeros(sc(1),num_sessions,sc(2),sc(3),sc(4),'single');
            stdcope=inf(sc(1),num_sessions,sc(2),sc(3),sc(4),'single');

        end;

        for con=1:size(first_level_results.cope,3),
            % check for infs or zeros in stdcope
            if sum(squash(isinf(first_level_results.stdcope(:,:,con,:))))>0,
                warning(['Session ' num2str(sesi) ' has some infinite stdcopes for first-level contrast ' num2str(con) '. This may be due to there being zero trials for one of the first level conditions.']);
            end;

            if any(squash(first_level_results.stdcope(:,:,con,:))==0),
                warning(['Session ' num2str(sesi) ' has some zero stdcopes for first-level contrast ' num2str(con) '. This may be due to there being zero trials for one of the first level conditions. Setting these to inf and continuing.']);

                XX=first_level_results.stdcope(:,:,con,:);
                XX(XX == 0) = inf;
                first_level_results.stdcope(:,:,con,:) = XX;
            end;
        end;

        cope(:,sesi,:,:,:)=first_level_results.cope; % nvox x nsub x ntpts x ncons x nfreqs

        stdcope(:,sesi,:,:,:)=first_level_results.stdcope;

        first_level_results= rmfield(first_level_results,'cope');
        first_level_results= rmfield(first_level_results,'stdcope');

        if(sesi==1),

            subject_level_results=first_level_results;
            subject_level_results.fname=[oat.first_level.name '_' oat.subject_level.name];
            subject_level_results.recon_method=first_level_results.recon_method;

            if ~strcmp(subject_level_results.recon_method,'none'),

                if(isfield(oat.first_level,'mni_coords'))

                    current_level_mni_coords         = oat.first_level.mni_coords;
                    subject_level_results.mni_coords = first_level_results.mni_coords;

                else

                    % use exact same coords/mask as first_level
                    subject_level_results.gridstep=first_level_results.gridstep;
                    S=[];

                    S.lower_level_mask_fname=[oat.source_recon.dirname '/' oat.first_level.name '_mask.nii.gz'];
                    S.current_level_mask_fname=[oat.source_recon.dirname '/' oat.first_level.name '_' oat.subject_level.name '_mask.nii.gz'];

                    S.current_level=oat.subject_level;
                    S.lower_level_mni_coord=first_level_results.mni_coords;
                    S.lower_level_gridstep=first_level_results.gridstep;

                    [subject_level_results.mask_indices_in_lower_level,subject_level_results.mni_coords]=setup_mask_indices(S);

                    clear S
                end;

            end;
        end;

    end;

    subject_level_results.subject_name=['subject' num2str(sub)];

    % voxel loop
    ft_progress('init', 'etf');

    if(size(cope,2)>1)
        if(isfield(oat.subject_level,'design_matrix')),

            for iVox = 1:size(cope,1),
                iVox
                ft_progress(iVox/size(cope,1));
                % time loop
                for iTime = 1:size(cope,3)
                    % contrast loop
                    for iContrast = 1:size(cope,4)
                        % freq loop
                        for iFreq = 1:size(cope,5)

                            cope_tmp=permute(cope(iVox,:,iTime,iContrast,iFreq),[2 1 3 4 5]);
                            stdcope_tmp_sq=permute(stdcope(iVox,:,iTime,iContrast,iFreq),[2 1 3 4 5]);

                            S=diag(stdcope_tmp_sq.^2);

                            z=oat.subject_level.design_matrix;
                            %z=ones(length(S),1);

                            [gam, covgam] = fixed_effects(cope_tmp,z,S);

                            cope_tmp=permute(cope,[2 1 3 4 5]);

                            subject_level_results.cope2(iVox,iTime,iContrast,iFreq)=gam;
                            subject_level_results.stdcope2(iVox,iTime,iContrast,iFreq)=sqrt(covgam);

                        end % for iFreq = 1:nFreqs
                    end % for iContrast = 1:nflips
                end % for iTime = 1:nTimes
            end % for iVox = 1:nVox

        else

            if any(isinf(stdcope))
                keyboard
            end

            cope_tmp=permute(cope,[2 1 3 4 5]);
            stdcope_tmp_sq=permute(stdcope,[2 1 3 4 5]).^2;

            recip_stdcope=1./stdcope_tmp_sq;
            covgam = 1./sum(recip_stdcope,1);
            temp = cope_tmp.* recip_stdcope;
            temp = sum(temp,1);
            gam = covgam.*temp;
            subject_level_results.cope = permute(gam,[2 3 4 5 1]);
            subject_level_results.stdcope = permute(sqrt(covgam),[2 3 4 5 1]);

        end;

    else

        disp('Only one session for subject - just passing stats through');
        subject_level_results.cope=permute(cope(:,1,:,:,:),[1 3 4 5 2]);
        subject_level_results.stdcope=permute(stdcope(:,1,:,:,:),[1 3 4 5 2]);

    end;

    ft_progress('close');

    %% Create contra-ipsi if requested
    % combine the COPE for one condition with a left-right flipped COPE from
    % another condition. This is useful for computing Ipsilateral Vs Contralateral
    % conditions by combining COPEs from stimulus left and stimulus right
    % conditions for instance.

    if ~isempty(oat.subject_level.merge_contraipsi)

        disp(['Doing fixed effects contra-ipsi merge for subject ' num2str(sub)]);

        
        [nflips, nconds] = size(oat.subject_level.merge_contraipsi);

        if nconds ~= 2
            error('Only two subject level copes can be merged!')
        end

        for iFlip = 1:nflips

            % Take a copy of the copes
            cope2    = subject_level_results.cope;
            stdcope2 = subject_level_results.stdcope;

            new_cope    = zeros(size(cope2,2),size(cope2,1),1,1,2);
            new_stdcope = zeros(size(cope2,2),size(cope2,1),1,1,2);

            cope_tmp    = permute(cope2,[2,1,3,4,5]);
            stdcope_tmp_sq = permute(stdcope2,[2,1,3,4,5]).^2;

            % Main loop
            ft_progress('init','etf');
            for iTime = 1:size(cope_tmp,1)
                ft_progress(iTime/size(cope_tmp,1));

                % Get contrast indices
                contrast_to_flip = oat.subject_level.merge_contraipsi(iFlip,1);
                contrast_not_to_flip = oat.subject_level.merge_contraipsi(iFlip,2);

                % Flip cope
                [flip_mni, flip_cope] = flip_voxels(first_level_results.mni_coords,cope_tmp(iTime,:,contrast_to_flip)');
                % Flip stdcope
                [flip_mni, flip_stdcope] = flip_voxels(first_level_results.mni_coords,stdcope_tmp_sq(iTime,:,contrast_to_flip)');

                % combine into new cope with a length 2 session
                % dimension on the end
                new_cope(iTime,:,:,:,:) = cat(5,flip_cope',cope_tmp(iTime,:,contrast_not_to_flip));
                new_stdcope(iTime,:,:,:,:) = cat(5,flip_stdcope',stdcope_tmp_sq(iTime,:,contrast_not_to_flip));
            end

            ft_progress('close');

            % perm to  [sessions, voxels, time,contrast, frequency]
            new_cope = permute(new_cope,[5,2,1,3,4]);
            new_stdcope = permute(new_stdcope,[5,2,1,3,4]);

            % fixed effects
            recip_stdcope=1./new_stdcope;
            covgam = 1./sum(recip_stdcope,1);
            temp = new_cope.* recip_stdcope;
            temp = sum(temp,1);
            gam = covgam.*temp;

            % Append our new contrast
            subject_level_results.cope(:,:,end+1) = permute(gam,[2 3 1]);
            subject_level_results.stdcope(:,:,end+1) = permute(sqrt(covgam),[2 3 1]);

            subject_level_results.first_level_contrast_name{end+1} = sprintf('ContraIpsi_cond:%s_cond:%s',num2str(contrast_to_flip),num2str(contrast_not_to_flip));

        end
    end

    %% Create laterality contrast COPE if requested
    % Take a fixed effects difference between a condition and a left-right
    % flipped version of itself. This is useful for computing differences between
    % responses in the two hemispheres. The resulting volume will be
    % symmetrical in the x axis (left-right) with values indicating the
    % difference between hemispheres

    if ~isempty(oat.subject_level.compute_laterality)

        disp(['Doing fixed effects laterality estimate for subject ' num2str(sub)]);
        
        [nflips] = length(oat.subject_level.compute_laterality);

        for iFlip = 1:nflips

            % Take a copy of the copes
            cope2    = subject_level_results.cope;
            stdcope2 = subject_level_results.stdcope;

            % [time, voxels, contrast, frequency, session]
            new_cope    = zeros(size(cope2,2),size(cope2,1),1,1,2);
            new_stdcope = zeros(size(cope2,2),size(cope2,1),1,1,2);

            cope_tmp    = permute(cope2,[2,1,3,4,5]);
            stdcope_tmp_sq = permute(stdcope2,[2,1,3,4,5]).^2;

            % Main loop
            ft_progress('init','etf');
            for iTime = 1:size(cope_tmp,1)
                ft_progress(iTime/size(cope_tmp,1));

                % Get contrast indices, these are the same here
                contrast_to_flip = oat.subject_level.compute_laterality(iFlip);
                contrast_not_to_flip = oat.subject_level.compute_laterality(iFlip);

                % Flip cope
                [flip_mni, flip_cope] = flip_voxels(first_level_results.mni_coords,cope_tmp(iTime,:,contrast_to_flip)');
                % Flip stdcope
                [flip_mni, flip_stdcope] = flip_voxels(first_level_results.mni_coords,stdcope_tmp_sq(iTime,:,contrast_to_flip)');

                % combine into new cope with a length 2 session
                % dimension on the end
                new_cope(iTime,:,:,:,:) = cat(5,flip_cope',cope_tmp(iTime,:,contrast_not_to_flip));
                new_stdcope(iTime,:,:,:,:) = cat(5,flip_stdcope',stdcope_tmp_sq(iTime,:,contrast_not_to_flip));
            end

            ft_progress('close');

            % perm to  [sessions, voxels, time,contrast, frequency]
            new_cope = permute(new_cope,[5,2,1,3,4]);
            new_stdcope = permute(new_stdcope,[5,2,1,3,4]);

            % fixed effects
            recip_stdcope=1./new_stdcope;
            covgam = 1./sum(recip_stdcope,1);
            temp = new_cope.* recip_stdcope;
            %temp = sum(temp,1);
            temp = temp(1,:,:) - temp(2,:,:);
            gam = covgam.*temp;

            % Append our new contrast
            subject_level_results.cope(:,:,end+1) = permute(gam,[2 3 1]);
            subject_level_results.stdcope(:,:,end+1) = permute(sqrt(covgam),[2 3 1]);

            subject_level_results.first_level_contrast_name{end+1} = sprintf('Laterality_cond:%s',num2str(contrast_to_flip));

        end
    end


    %% save results to file
    subject_level_results.S=oat.subject_level;
    subject_level_results.first_level_name=first_level_results.name;

    subject_level_results.subject_level=oat.subject_level;

    subject_level_results.name=oat.subject_level.name;
    subject_level_results.recon_method=first_level_results.recon_method;

    subject_level_results.fname=[ subject_level_results.subject_name '_' oat.first_level.name '_' oat.subject_level.name ];
    disp(['Saving statistics in file ' oat.source_recon.dirname '/' subject_level_results.fname]);
    oat_save_results( oat, subject_level_results );
    results_fnames{sub}=subject_level_results.fname;

end;
