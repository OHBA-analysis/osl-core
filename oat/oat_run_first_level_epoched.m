function [ results_fnames first_level_results ] = oat_run_first_level_epoched( oat )

% [ results_fnames ] = oat_run_first_level_epoched( oat )
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct
% and runs first level trialwise GLM
%
% This function should normally be called using osl_run_oat(oat);
%
% MWW 2012

OSLDIR = getenv('OSLDIR');

do_glm_statewise=oat.first_level.hmm_do_glm_statewise;

first_level=oat.first_level;
source_recon_name='source_recon';

% check if user has asked not to do GLM
if(~isfield(first_level,'doGLM')),
    first_level.doGLM = 1;
end

if(~isfield(first_level,'bc_trialwise')),
    first_level.bc_trialwise=0;
end;

if(~isfield(first_level, 'trialwise_directory')),
    first_level.trialwise_directory = [];
end

if isfield(first_level,'time_downsample_factor') && ~isempty(first_level,'time_downsample_factor'),
    error('first_level.time_downsample_factor no longer supported');
end;

if isfield(first_level,'tf_downsample_factor') && ~isempty(first_level,'tf_downsample_factor'),
    error('first_level.tf_downsample_factor no longer supported');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the downsample factors
% tf_downsample_factor is used in source space after recon weights are applied
% time_downsample_factor is used after TF transform (as part of osl_tf_transform call) in sensor space on complex values
if ~isempty(first_level.post_tf_downsample_factor)
    post_tf_ds_factor = 1./first_level.post_tf_downsample_factor;
else
    post_tf_ds_factor = 1;
end

if ~isempty(first_level.post_movingaverage_downsample_factor)
    post_movingaverage_ds_factor = 1./first_level.post_movingaverage_downsample_factor;
else
    post_movingaverage_ds_factor = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set contrasts

contrast_list=first_level.contrast;
for c=1:length(contrast_list),
    contrast_list{c}=contrast_list{c}(:);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check some of the settings

if(sum(first_level.bc)>0 && first_level.baseline_timespan(1)>first_level.baseline_timespan(2)),
    error('first_level.baseline_timespan(1)>first_level.baseline_timespan(2)');
end;

if(~isfield(oat.first_level,'design_matrix') && ~isfield(oat.first_level,'design_matrix_summary')),
    error('Design matrix is not specified');
end;

if(length(oat.first_level.bc)~=length(oat.first_level.contrast)),
    error('first_level.bc and first_level.contrasts need to be the same length');
end;

if(length(oat.first_level.contrast_name)~=length(oat.first_level.contrast))
    error('oat.first_level.contrast and oat.first_level.contrast_name need to be same length');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set first level diagnostic report up
report_dir=[oat.results.plotsdir '/' oat.first_level.name];
first_level_results.report=osl_report_setup(report_dir,'First level (epoched)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over sessions
clear stats;
for subi_todo=1:length(first_level.sessions_to_do),

    sub=first_level.sessions_to_do(subi_todo);

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OAT FIRST LEVEL ON SESS = ' num2str(sub) '  %%%%%%%%%%%%%%%%%%%%%%%'])

    %sub/length(first_level.submatfiles),

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load in the beamformer result for this session
    source_recon_results_fname = oat.source_recon.results_fnames{sub};
    source_recon_results=oat_load_results(oat,source_recon_results_fname);

    % get D and update path in case OAT dir has been moved
    D = oat_get_sensordata( source_recon_results );
    D.fullfile
    % Make sure we don't stamp on the source recon
    S = [];
    S.D = D;
    tmp = char(D.fullfile);
    S.outfile = [tmp(1:end-4) '_firstlevel'];
    D = spm_eeg_copy(S);
    D.fullfile
    %% setup things that are common to all sessions

    %if(subi_todo==1),
    % This loop has been commented out as some later operations change array
    % sizes further down the analysis path, for instance
    % first_level_results.cope has 204 channels rather than 305 after grads
    % have been combined. This then breaks subsequent iterations which will
    % have the orginial number of channels. There is probably a better fix. AQ

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% set up mask indices from mask/coords
        if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis

            disp('Working in sensor space');

            is_sensor_space=1;

            if strcmp(oat.source_recon.modalities{1},'EEG')   % added by DM
                modality_meeg='EEG';
            else
                modality_meeg='MEG';
            end
            chanindmeg = strmatch(modality_meeg, D.chantype);

            first_level_results.chanind=chanindmeg;
            first_level_results.chanlabels=D.chanlabels(chanindmeg);
            first_level_results.chantype=D.chantype(chanindmeg);

            first_level_results.mask_indices_in_source_recon=chanindmeg;

        else,

            is_sensor_space=0;

            if(isfield(first_level,'mni_coords'))

                current_level_mni_coord = first_level.mni_coords;

                % load in lower level coords
                lower_level_mni_coord=source_recon_results.mni_coords;
                %lower_level_mask_fname  = [oat.source_recon.dirname '/' source_recon_name '_mask'];
                %[ lower_level_mni_coord, ~ ] = osl_mnimask2mnicoords(lower_level_mask_fname);

                indices_in_lower_level=zeros(size(current_level_mni_coord,1),1);
                for vox=1:size(current_level_mni_coord,1),
                    [indices_in_lower_level(vox), vec, dist]=nearest_vec(lower_level_mni_coord,current_level_mni_coord(vox,:));

                    if(dist>source_recon_results.gridstep),
                        warning(['Distance between desired coordinate [' num2str(current_level_mni_coord(vox,:)) '] and nearest in source_recon results is ' num2str(dist) 'mm.']);
                    end;
                end;

                first_level_results.mask_indices_in_source_recon = indices_in_lower_level;
                first_level_results.mni_coords = lower_level_mni_coord(indices_in_lower_level,:);

                % the mni coordinates from the specified list
                str=['Desired MNI coordinates ' ];
                for vox=1:size(current_level_mni_coord,1),
                    str=[str ', [' num2str(current_level_mni_coord(vox,:)) ']'];
                end;
                disp(str);

                % the nearest mni coordinates from the lower level mask
                str=['Nearest dipoles from lower level at MNI coordinates ' ];
                for vox=1:length(indices_in_lower_level),
                    str=[str ', [' num2str(lower_level_mni_coord(indices_in_lower_level(vox),:)) ']'];
                end;
                disp(str);

            elseif ~first_level.parcellation.do
                % setup std space brain
                first_level_results.gridstep=source_recon_results.gridstep;

                S=[];
                S.lower_level_mask_fname=[oat.source_recon.dirname '/' source_recon_name '_mask.nii.gz'];
                S.current_level_mask_fname=[oat.source_recon.dirname '/' oat.first_level.name '_mask.nii.gz'];
                S.current_level=first_level;
                S.lower_level_mni_coord=source_recon_results.mni_coords;
                S.lower_level_gridstep=source_recon_results.gridstep;

                [first_level_results.mask_indices_in_source_recon, first_level_results.mni_coords]=setup_mask_indices(S);
                lower_level_mni_coord=source_recon_results.mni_coords;

                clear S;
            else
                if first_level.parcellation.do

                    first_level_results.gridstep=source_recon_results.gridstep;

                    Dold=D;
                    S=first_level.parcellation;
                    S.prefix='p';
                    S.D=D;
                    D = osl_apply_parcellation(S);

                    %%%%%%%%%%%%%%%%%%%
                    % add back in class channel
                    classchanind=find(strcmp(Dold.chanlabels,'Class'));

                    Sc=[];
                    Sc.D=D;
                    Sc.newchandata=[Dold(classchanind,:,:)];
                    Sc.newchanlabels{1}='Class';
                    Sc.newchantype{1}='CLASS';

                    [ Dnew ] = osl_concat_spm_eeg_chans( Sc );

                    D.delete;
                    D=Dnew;
                    %%%%%%%%%%%%%%%%%%%

                    first_level_results.mask_indices_in_source_recon=1:size(D,1)-1;
                    first_level_results.mni_coords=D.parcellation.mni_coords;
                end

            end;

        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% set up time and frequency info by making dummy call to osl_tf_transform
        % Note that any time or frequency averaging will be done AFTER the time-frequency decomposition

        if(~isempty(first_level.tf_hilbert_freq_ranges)),
            disp('oat.first_level.tf_hilbert_freq_ranges is set, freq_range, tf_freq_res and tf_num_freqs will be overridden');

            first_level.tf_num_freqs=size(first_level.tf_hilbert_freq_ranges,1);
            freq_range=[];
        else
            freq_range = first_level.tf_freq_range; % range used for stat testing
            if isempty(freq_range)
                % if not set as a first_level option, default to the range used
                % for source_recon
                freq_range=source_recon_results.source_recon.freq_range;

                if isempty(freq_range)
                    freq_range=[1 D.fsample/3];
                    disp(['WARNING: No frequency range has been specified. Set oat.source_recon.freq_range or oat.first_level.tf_freq_range. Using ' num2str(freq_range) ' Hz.']);
                end % isempty(freq_range)
            end % if isempty(freq_range)

        end;

        if ~isempty(freq_range),
            if(freq_range(2)>D.fsample/2),
                warning(['Nyquist=' num2str(D.fsample/2) ', freq_range(2)>D.fsample/2. Results will be dodgey due to looking at frequencies higher than the Nyquist limit.']);
            end;
        else
            if(max(squash(first_level.tf_hilbert_freq_ranges))>D.fsample/2),
                warning(['Nyquist=' num2str(D.fsample/2) ', freq_range(2)>D.fsample/2. Results will be dodgey due to looking at frequencies higher than the Nyquist limit.']);
            end;
        end;


        D_time_indices=find(source_recon_results.samples2use); % time indices in D.time used by source recon
        D_times = D.time(D_time_indices); % times for time window of data that will be the input into the glm

        if isempty(first_level.time_range)
            first_level.time_range=source_recon_results.woi;
        end

        tf_settings = [];
        tf_settings.tf_method               = first_level.tf_method;
        tf_settings.tf_hilbert_freq_res     = first_level.tf_hilbert_freq_res;
        tf_settings.tf_hilbert_freq_ranges  = first_level.tf_hilbert_freq_ranges;
        tf_settings.tf_freq_range           = freq_range; % set above
        tf_settings.tf_num_freqs            = first_level.tf_num_freqs;
        tf_settings.tf_hanning_ncycles      = first_level.tf_hanning_ncycles;
        tf_settings.raw_times               = D_times;
        tf_settings.ds_factor               = post_tf_ds_factor;
        tf_settings.tf_morlet_factor        = first_level.tf_morlet_factor;
        tf_settings.tf_calc_amplitude=0; % we will calc the amplitude after applying recon weights
        tf_settings.tf_hilbert_do_bandpass_for_single_freq = first_level.tf_hilbert_do_bandpass_for_single_freq;

        out = osl_tf_transform( tf_settings ); % passed without a data matrix, osl_tf_transform just gets container sizes

        if isfield(out, 'tf_morlet_basis')
            first_level.tf_morlet_basis = out.tf_morlet_basis;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% setup time and frequency windows/indices for using later

        % setup result time dimension
        %tf_time_indices_into_D_times = ( out.tf_times >= first_level.time_range(1) ) & ( out.tf_times <= first_level.time_range(2) );
        tf_time_indices_into_D_times = false(1,numel(out.tf_times));
        for i=1:size(first_level.time_range,1)
            tf_time_indices_into_D_times(( out.tf_times >= first_level.time_range(i,1) ) & ( out.tf_times <= first_level.time_range(i,2) )) = true;
        end

        first_level_results.times = out.tf_times(tf_time_indices_into_D_times);

        tf_out_times=out.tf_times(tf_time_indices_into_D_times);

        % downsample first_level_results.times
        % temporally downsample
        if post_movingaverage_ds_factor~=1,
            tmp=ones(size(first_level_results.times));
            tmp=osl_spm_resample(tmp,post_movingaverage_ds_factor);
            first_level_results.times=linspace(first_level_results.times(1),first_level_results.times(end),length(tmp));
        end;

        if first_level.time_average
            first_level_results.times=mean(first_level_results.times);
            tres=1;
        end

        % setup result freq dimension
        if ~isfield(out,'tf_freq_ranges')
            out.tf_freq_ranges=repmat(out.tf_freqs,2,1)';
            if isfield(out,'tf_freq_res')
                out.tf_freq_ranges(:,1) = out.tf_freq_ranges(:,1) - out.tf_freq_res/2;
                out.tf_freq_ranges(:,2) = out.tf_freq_ranges(:,2) + out.tf_freq_res/2;
            end;
        end

        first_level_results.frequency_ranges=out.tf_freq_ranges;
        first_level_results.frequencies = out.tf_freqs;

        % setup result spatial dimension
        Nvoxels_out=length(first_level_results.mask_indices_in_source_recon);

        clear out

        %%%%%%%%%%%%%%%%
        %% results containers
        ntpts                       = length(first_level_results.times);
        nfreqs                      = length(first_level_results.frequencies);

        artefact_map = nan(Nvoxels_out,nfreqs);

        if(first_level.doGLM)
            first_level_results.stdcope=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs]);% num_voxels x  num_timepoints x num_contrasts x num_freqs
            first_level_results.cope=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs]);% num_voxels x num_timepoints x num_contrasts x num_freqs

        end;

        first_level_results.pseudo_zstat_var = zeros([Nvoxels_out,1],'single');

    %end; % end of if(sub==1)

    source_recon_results.is_sensor_space=is_sensor_space;

    disp(['Reconstruct time courses and computing stats for dataset ' source_recon_results.fname]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% setup HMM
    if strcmp(source_recon_results.recon_method,'none'),
        NKrecon=1;
        NKglm=NKrecon;
    else
        NKrecon=1;%numel(source_recon_results.BF.inverse.MEG.class);
        NKglm=NKrecon;

        if NKrecon>1 && ~do_glm_statewise
            warning('Must do GLM statewise when using multiple classes (e.g. HMM)');
        end;

        if(~do_glm_statewise),
            NKglm=1;
        end;

        if(NKglm>1)
            if post_movingaverage_ds_factor~=1,
                error('Can not do any downsampling when using multiple classes in GLM (e.g. HMM)');
            end;
        end;
    end;

    if(first_level.doGLM && NKglm>1)
        first_level_results.stdcope_by_state=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs,NKglm]);% num_voxels x  num_timepoints x num_contrasts x num_freqs x num_states
        first_level_results.cope_by_state=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs,NKglm]);% num_voxels x num_timepoints x num_contrasts x num_freqs x num_states
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% setup the GLM design matrix for this session
    % setup trial types
    trialtypes=[];
    triallist=[];

    for i=1:length(source_recon_results.source_recon.conditions), % indexes conditions/triggers within subset

        trigname=source_recon_results.source_recon.conditions{i};
        Ntrialspercond=length(D.indtrial(trigname,'GOOD')); %% number of trials for this condition
        triallist=[triallist , D.indtrial(trigname,'GOOD')];
        trialtypes=[trialtypes ; ones(Ntrialspercond,1).*i];

    end;

    S3=[];
    S3.trialtypes=trialtypes;

    if numel(first_level.design_matrix_summary)==length(oat.source_recon.results_fnames) &&(ischar(first_level.design_matrix_summary{sub})),

        S3.Xsummary=first_level.design_matrix_summary{sub};

        if(isfield(first_level,'trial_rejects'))
            S3.trial_rejects=first_level.trial_rejects;
        else
            warning('oat.first_level.trial_rejects not set');
        end;
    else
        S3.Xsummary=first_level.design_matrix_summary;
    end;

    first_level.x=oat_setup_designmatrix(S3);

    const_reg=0;
    for p=1:size(first_level.x,2),
        if(range(first_level.x(:,p))==0),
            if(const_reg),
                warning('Multiple constant regressors in the design matrix!!!!?');
            end;
            const_reg=1;

        end;

        if(~any(first_level.x(:,p))),
            warning(['Empty (all zeros) regressor (EV' num2str(p) ') in the design matrix!!!!? Could cause problems for subsequent OLS group analysis.']);
        end;
    end;

    if(first_level.do_glm_demean)
        disp('CAREFUL: are you sure you want the first_level.do_glm_demean flag off? This means that main effect contrasts will not be interpetable as mean effects');
    end;

    if(const_reg && first_level.do_glm_demean)
        disp('CAREFUL: are you sure you want the first_level.do_glm_demean flag on with constant regressors in the design matrix!!!!?');
    elseif(~const_reg && ~first_level.do_glm_demean)
        disp('CAREFUL: are you sure you want no first_level.do_glm_demean flag on with no constant regressors in the design matrix!!!!?');
    end;

    Ntrials   = length(triallist);

    %%%%%%%%%%%%%%%
    %% load in sensor data
    first_level_results.D_sensor_data=D;

    classchanind=find(strcmp(D.chanlabels,'Class'));
    if(isempty(classchanind)),
        error(['No ''Class'' chanlabel in: ' D.fname]);
    end;

    % class_samples_inds_recon is used for using weights in recon
    for si=1:NKrecon,
        class_samples_inds_recon{si} = (D(classchanind, D_time_indices, triallist)==si);
    end;

    % class_samples_inds_glm is used for doing GLM state(class)-wise
    if(~do_glm_statewise),
        class_samples_inds_glm{1}=ones(1,length(D_time_indices),length(triallist));
    else
        class_samples_inds_glm=class_samples_inds_recon;
    end;

    for si=1:NKrecon,
        state_tinds=class_samples_inds_recon{si};

        total_time=(D_times(end)-D_times(1))*size(state_tinds,3);
        state_time=total_time*sum(sum(squeeze(state_tinds)))/prod(size(squeeze(state_tinds)));
        txt=['State ' num2str(si) ' is active for ' num2str(state_time) 'secs'];

        if(state_time<10),
            warning([txt '. It will be excluded from the GLM analysis.']);
        else
            disp(txt);
        end;
    end;
    clear state_tinds;

    %%%%%%%%%%%%%%%
    %% setup results containers
    if(first_level.doGLM)
        disp('First level COPEs outputted will have dimension Nvoxels x Ntpts x Ncontrasts x Nfreqs:');
        tmp=[size(first_level_results.cope) 1 1];
        disp([num2str(tmp(1:4))]);
    end;

    % containers
    if(first_level.doGLM)
        covarpe_by_state=zeros([ntpts,size(first_level.x,2),size(first_level.x,2),length(first_level_results.frequencies),NKglm]);% num_timepoints x num_contrasts x num_contrasts x num_freqs x num_states
        pe_by_state=zeros([ntpts,size(first_level.x,2),length(first_level_results.frequencies),NKglm]);% num_timepoints x num_contrasts x num_freqs x num_states
    end;

    if first_level.save_trialwise_data

        if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis
            trlwise.labels = first_level_results.chanlabels;
        else
            for ilab = 1:Nvoxels_out; trlwise.labels{ilab} = ['voxel' num2str(ilab)]; end;
        end % strcmp(source_recon_results.recon_method,'none'), % sensor space analysis

            % transform-specific parameters
            switch first_level.tf_method
                case 'none'
                    trlwise.type    = 'time_domain';
                    trlwise.dat     = nan(Nvoxels_out,Ntrials,ntpts,nfreqs);
                    trlwise.times   = first_level_results.times;
                otherwise
                    trlwise.type    = first_level.tf_method;
                    trlwise.dat     = nan(Nvoxels_out,Ntrials,ntpts,nfreqs);
                    trlwise.freqs   = first_level_results.frequencies;
                    trlwise.times   = first_level_results.times;
            end
            trlwise.triallist = triallist;
            trlwise.rawdat    = fullfile(D.path,D.fname);
            %trlwise.nonblinktimes = class_samples_inds_glm_tf{1};
            trlwise.fsample = 1/mode(diff(D_times));

            % spatial parameters
            isSensorSpace = (isfield(source_recon_results, 'is_sensor_space') && ...
                             ~source_recon_results.is_sensor_space)           || ...
                            strcmpi(source_recon_results.recon_method, 'none');
            if ~isSensorSpace,
                if (isfield(first_level,'mni_coords'))
                    trlwise.mni_coords   = first_level_results.mni_coords;
                else
                    trlwise.mask         = nii.load([oat.source_recon.dirname '/' oat.first_level.name '_mask']);
                    trlwise.gridstep     = first_level_results.gridstep;
                    trlwise.mni_coords   = first_level_results.mni_coords;
                end % if (isfield(first_level,'mni_coords'))
            end;

    end % first_level.trialwise_data

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% do outer loop over freqs, this helps limit RAM usage when doing TF transforms
    for f=1:first_level.tf_num_freqs,

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Do T-F transform in sensor space
        %% Transform the data - either downsampling in time domain, or a TF transform

        if(~strcmp(first_level.tf_method,'none'))

            disp(['Doing T-F transform in sensor space on ' mat2str(first_level_results.frequency_ranges(f,:))  ' Hz...']);
        end;

        tf_settings = [];
        tf_settings.tf_method               = first_level.tf_method;
        tf_settings.tf_hilbert_freq_res     = first_level.tf_hilbert_freq_res;
        tf_settings.tf_hilbert_freq_ranges  = first_level_results.frequency_ranges(f,:);
        tf_settings.tf_freq_range           = first_level_results.frequency_ranges(f,:); % set above
        tf_settings.tf_num_freqs            = 1;
        tf_settings.raw_times               = D_times;
        tf_settings.ds_factor               = post_tf_ds_factor;
        tf_settings.tf_morlet_factor        = first_level.tf_morlet_factor;
        tf_settings.tf_hanning_ncycles      = first_level.tf_hanning_ncycles;
        if is_sensor_space == 1
            tf_settings.tf_calc_amplitude = 1;
        else
            % This is broken, the spm-object conversion below (D_tf =
            % osl_change_spm_eeg_data( Sc ); % line 638) drops the imaginary part of
            % the TF information. This is as clone sliently preserves the
            % input data-type and doesn't give a warning.
            %
            % Going to change this to output the amplitude, if we want to
            % do this afte the weights we will need a workaround for the
            % data copy and to uncomment the amplitude computation on line
            % 828
            %tf_settings.tf_calc_amplitude=0; % we will calc the amplitude after applying recon weights
            tf_settings.tf_calc_amplitude=1; % we will calc the amplitude now
        end
        if nfreqs>1,
            tf_settings.tf_hilbert_do_bandpass_for_single_freq = 1;
        else
            tf_settings.tf_hilbert_do_bandpass_for_single_freq = first_level.tf_hilbert_do_bandpass_for_single_freq;
        end;

        if isfield(first_level,'tf_morlet_basis')
            tf_settings.tf_morlet_basis={};
            tf_settings.tf_morlet_basis{1} = first_level.tf_morlet_basis{f};
        end

        sensor_data=D(:, D_time_indices, triallist);

        dims=size(sensor_data);
        dims_tf=[dims(1),length(tf_time_indices_into_D_times),Ntrials,tf_settings.tf_num_freqs];
        disp(['size(sens_data_tf)=' mat2str((dims_tf)) ', (nsens x ntpts x ntri x nfreq)']);

        if(prod(size(dims_tf)) > 1e9)
            warning(['size(sens_data_tf)=' mat2str((dims_tf)) '. May cause RAM problems. Consider using less freq bins or using oat.first_level.post_tf_downsample_factor>1.']);
        end;

        out = osl_tf_transform( tf_settings, sensor_data);
        sensor_data_tf = out.dattf; % [channel, time, trial, frequency]

        % get indices into the class vectors, in order to downsample as per
        % osl_tf_transform
        ds_inds = round(linspace(1,numel(D_times),numel(out.tf_times)));
        for si=1:NKglm,
            class_samples_inds_glm_tf{si}   = class_samples_inds_glm{si}(:,ds_inds,:);
        end;
        for si=1:NKrecon,
            class_samples_inds_recon_tf{si} = class_samples_inds_recon{si}(:,ds_inds,:);
        end;

        %if(first_level.tf_num_freqs>1), ft_progress('close');end;

        %%%%%%%%%%%%%%%%
        %% Now do windowing over time
        sensor_data_tf=sensor_data_tf(:,tf_time_indices_into_D_times,:,:);
        for si=1:NKglm,
            class_samples_inds_glm_tf{si}=class_samples_inds_glm_tf{si}(:,tf_time_indices_into_D_times,:);
        end;
        for si=1:NKrecon,
            class_samples_inds_recon_tf{si}=class_samples_inds_recon_tf{si}(:,tf_time_indices_into_D_times,:);
        end;

        %%%%%%%%%%%%%%%%
        %% create new D object from sensor_data_tf but which is otherwise a
        % copy of original D object
        % i.e. following any TF transform and time windowing
        % note that this will ony contain one freq band as we are looping
        % over freq bands
        Sc=[];
        Sc.D = D;
        [path nm ext]=fileparts(fullfile(D));
        Sc.newname = [path '/TF' nm ext];
        Sc.newdata = sensor_data_tf;
        Sc.time = tf_out_times;
        Sc.frequencies = first_level_results.frequencies(f);
        Sc.remove_montages=0;
        Sc.cond_list = D.conditions(triallist);
        D_tf = osl_change_spm_eeg_data( Sc );

        % add back in Class channel:
        classchanind=find(strcmp(D.chanlabels,'Class'));
        classchanind_tf=find(strcmp(D_tf.chanlabels,'Class'));
        D_tf(classchanind_tf,:,:)=D(classchanind,D_time_indices(tf_time_indices_into_D_times),triallist);

        clear sensor_data_tf;

        %%%%%%%%%%%%%%%%
        %% Temporally downsample state time courses (data will be done later after taking abs (for t-f data) inside the voxel loop)
        for si=1:NKglm,

            state_tinds=class_samples_inds_glm_tf{si};

            % moving average within the window specified
            if ~isempty(first_level.time_moving_av_win_size)

                if(NKglm>1)
                    error('first_level.time_moving_av_win_size not compatible with more than one class');
                end;

                % do nothing
            end

            % temporally downsample
            if post_movingaverage_ds_factor~=1,
                clear state_tinds_new;
                for tri=1:size(state_tinds,3),
                    tmp=permute(state_tinds(1,:,tri),[1,2,3]);
                    tmp=osl_spm_resample(tmp,post_movingaverage_ds_factor);
                    state_tinds_new(1,:,tri)=tmp;
                end;
                state_tinds=state_tinds_new;

                % need binary states
                state_tinds=round(state_tinds);
                state_tinds(state_tinds<0)=0;
                state_tinds(state_tinds>1)=1;
            end;

            class_samples_inds_glm_tf{si}=state_tinds;
        end;

        % no need to temporally downsample or moving average class_samples_inds_recon_tf

        % downsample first_level_results.times


        %%%%%%%%%%%%%%%%
        %% setup design matrices that can be reused over voxels
        ntpts=length(first_level_results.times);
        x_stored=cell(NKglm,ntpts);
        %pinvxtx_stored=cell(NKglm,ntpts);
        %pinvx_stored=cell(NKglm,ntpts);
        contrast_list_mat=cell2mat(contrast_list);
        contrast_list_mat2=eye(size(first_level.x,2));

        isrepresented = zeros(NKglm,size(state_tinds,2));
        state_times=zeros(NKglm,1);

        for si=1:NKglm,

            state_tinds=class_samples_inds_glm_tf{si};

            total_time=(D_times(end)-D_times(1))*size(state_tinds,3);
            state_time=total_time*sum(sum(squeeze(state_tinds)))/prod(size(squeeze(state_tinds)));
            state_times(si)=state_time;

            min_state_time=10;
            if(state_time>min_state_time),
                xc=[];
                contrast_list_matc=[];
                contrast_list_matc2=[];
                pinvxctxc_stored2{si}=[];
                pinvxc_stored2{si}=[];

                for t = 1 : size(state_tinds,2), % indexes time points within trials

                    %t=nearest(D_times,first_level_results.times(ts));

                    %disp(num2str(t/ntpts));

                    % x is (ntrials x nevs)
                    x=first_level.x;

                    % check if the state is associated with any data for this
                    % timepoint
                    if isempty(find(permute(state_tinds(1,t,:),[3 1 2])))
                    else
                        isrepresented(si,t) = 1;
                    end

                    x=x(find(permute(state_tinds(1,t,:),[3 1 2])),:);

                    if(oat.first_level.do_glm_demean)
                        for pp=1:size(x,2),
                            x(:,pp)=x(:,pp)-mean(x(:,pp));
                        end;
                    end;

                    x_stored{si,t}=x;

                    %pinvxtx_stored{si,t}=pinv(x'*x);
                    %pinvx_stored{si,t}=pinv(x);
                    pinvxtx_stored=pinv(x'*x);
                    pinvx_stored=pinv(x);

                    xc=sparse(blkdiag(xc,x));
                    pinvxctxc_stored2{si}=sparse(blkdiag(pinvxctxc_stored2{si},pinvxtx_stored));
                    pinvxc_stored2{si}=sparse(blkdiag(pinvxc_stored2{si},pinvx_stored));

                    contrast_list_matc=sparse(blkdiag(contrast_list_matc,contrast_list_mat));

                    contrast_list_matc2=sparse(blkdiag(contrast_list_matc2,contrast_list_mat2));

                end;


                contrast_list_matc=sparse(contrast_list_matc);
                contrast_list_matc2=sparse(contrast_list_matc2);

                % xc is (blkdiag of x repeated for each time point)
                % xc is ([ntpts*ntrials] x [ntpts*nevs])

                xc_stored{si}=sparse(xc);

                if(0)
                    xcxc=full(xc'*xc);
                    pinvxctxc_stored{si}=sparse(pinv(xcxc));

                    xc=full(xc);
                    pinvxc_stored{si}=sparse(pinv(xc));

                else
                    pinvxctxc_stored{si}=pinvxctxc_stored2{si};
                    pinvxc_stored{si}=pinvxc_stored2{si};
                end;
            end;
        end;

        if ~any(state_times>min_state_time),
            if NKglm==1,
                error(['Insufficient time points to proceed (<' num2str(min_state_time) 'secs)']);
            else
                error(['No states with sufficient time points to proceed (>' num2str(min_state_time) 'secs)']);
            end;
        end;

        clear state_tinds;

        % check that the non-artefact states account for all timepoints
        if ~isempty(oat.source_recon.artefact_chanlabel)
            iscaptured = sum(isrepresented(1:size(isrepresented,1)-1,:),1);
            if any(~iscaptured)
                percentmissing = 100 * (sum(~iscaptured) / numel(iscaptured));
                warning(['The non-artefact state(s) fail to cover ' num2str(percentmissing) '% of timepoints within the trials']);
            end
        end

        %%%%%%%%%%%%%%%%
        %% do everything from now on one voxel at a time (for the sake of
        %% RAM) in the first level mask space
        ft_progress('init', 'etf');

        % debug:
        %[indfind, mni_coord]=nearest_vec(source_recon_results.mni_coords,[-32 -26 50]);

        S2=[];
        S2.D=D_tf;
        if source_recon_results.is_sensor_space
            S2.montage_index=0;
        else
            if oat.first_level.do_weights_normalisation
                S2.montage_index=NKrecon+1; % with weights normalisation
            else
                S2.montage_index=1; % without weights normalisation
            end;
        end;
        S2.D_block.size=min(500,Nvoxels_out);
        first_level_results.sensor_D_tf=D_tf;

        for indind=1:Nvoxels_out, % indexes brain space

            ft_progress(indind/Nvoxels_out);

            S2.index=first_level_results.mask_indices_in_source_recon(indind);

            [dat_tf S2] = osl_get_recon_timecourse( S2 );

            % Place-holder amplitude computation - see lines 566
%             if tf_settings.tf_calc_amplitude == 0
%                 % Calculate power now
%                 dat_tf = (sqrt(dat_tf.*conj(dat_tf)));
%             end

            % dat_tf needs to trials x time:
            dat_tf=permute(dat_tf,[3 2 1]);

            first_level_results.pseudo_zstat_var(indind)=var(squash(dat_tf));

            if sum(isnan(squash(dat_tf)))~=0,

                if(first_level.doGLM)
                    first_level_results.cope(indind,:,:,:)=0;
                    first_level_results.stdcope(indind,:,:,:)=inf;
                    if(NKglm>1)
                        first_level_results.cope_by_state(indind,:,:,:,:)=0;
                        first_level_results.stdcope_by_state(indind,:,:,:,:)=inf;
                    end;
                end;
                error('NANs in recon data');

            else

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% calc the amplitude (up to now dat_tf is complex)
                if ~strcmp(first_level.tf_method,'none'),
                    dat_tf=abs(dat_tf);
                end;

                dat=dat_tf;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% moving average within
                %% the window specified
                if ~isempty(first_level.time_moving_av_win_size)

                    fsample=1/(D_times(2)-D_times(1));
                    windowSize = first_level.time_moving_av_win_size*fsample;

                    clear dat_new;
                    for tri=1:size(dat,1),
                            tmp=permute(dat(tri,:,1),[1,2,3]);
                            tmp=moving(tmp,round(windowSize));
                            dat_new(tri,:,1)=tmp;
                    end;
                    dat=dat_new;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% temporally downsample (needs to be done after abs (for t-f data))
                if post_movingaverage_ds_factor~=1,
                    clear dat_new;
                    for tri=1:size(dat,1),
                            tmp=permute(dat(tri,:,1),[1,2,3]);
                            tmp=osl_spm_resample(tmp,post_movingaverage_ds_factor);
                            dat_new(tri,:,1)=tmp;
                    end;
                    dat=dat_new;

                end;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% now do averaging over time
                if(first_level.time_average)
                    dat=mean(dat,2);
                end;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if first_level.save_trialwise_data
                    trlwise.dat(indind,:,:,f) = dat;
                end % if first_level.save_data

                if(first_level.bc_trialwise),
                    baseline_time_indices=first_level_results.times<first_level.baseline_timespan(2) & first_level_results.times>first_level.baseline_timespan(1);

                    if(sum(baseline_time_indices)>0)
                        % normalise by noise in baseline window:
                        for tri=1:size(dat,1), % indexes trials
                                dat(tri,:,1)=(dat(tri,:,1)-mean(dat(tri,baseline_time_indices,f)));
                        end;
                    else
                        warning('Could not baseline correct as there are no timepoints in the specified baseline period');
                    end;
                end;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Do GLM
                if first_level.doGLM,
                    npes=size(contrast_list_mat2,2);

                    if(size(dat,1)>1)
                            clear total_time state_time;
                            for si=1:NKglm,

                                state_tinds=class_samples_inds_glm_tf{si};

                                total_time{si}=(D_times(end)-D_times(1))*size(state_tinds,3);
                                state_time{si}=total_time{si}*sum(sum(squeeze(state_tinds)))/prod(size(squeeze(state_tinds)));

                                if(state_time{si}<10),

                                    if(NKglm>1)
                                        %first_level_results.cope_by_state(indind,:,:,f,si)=0;
                                        %first_level_results.stdcope_by_state(indind,:,:,f,si)=inf;

                                        pe_by_state(:,:,f,si)=0;
                                        for t= 1 : size(dat,2),
                                            covarpe_by_state(t,:,:,f,si)=inf;
                                        end;

                                    end;

                                else

                                    dattfc=zeros(size(dat,2)*size(state_tinds,3),1);
                                    %dattfc=[];
                                    iter=1;
                                    for t = 1 : size(dat,2), % indexes time points within trials

                                        dattf=dat(find(permute(state_tinds(1,t,:),[3 1 2])),t,1);

                                        if(oat.first_level.do_glm_demean)
                                            dattf=dattf-mean(dattf);
                                        end;

                                        %dattfc=cat(1,dattfc,dattf);
                                        dattfc(iter:iter+length(dattf)-1)=dattf;
                                        iter=iter+length(dattf);

                                    end;
                                    dattfc=dattfc(1:iter-1);

                                    % xc is (blkdiag of x from repeating x for each time point)
                                    % x is ([ntpts*ntrials] x [ntpts*nevs])
                                    x=xc_stored{si};
                                    pinvxtx=pinvxctxc_stored{si};
                                    pinvx=pinvxc_stored{si};

                                    [peout, covarpeout]=glm_fast_for_meg(dattfc,x,pinvxtx,pinvx);

                                    pe_by_state(:,:,f,si)=reshape(peout,npes,size(dat,2))';

                                    for t= 1 : size(dat,2),
                                        covarpe_by_state(t,:,:,f,si)=covarpeout((t-1)*npes+1:t*npes,(t-1)*npes+1:t*npes);
                                    end;

                                end;
                            end; % si (state index)

                            %% do fixed effects averaging of states
                            % time loop

                            if(NKglm>1)

                                % one state may be an artefact state
                                if ~isempty(oat.source_recon.artefact_chanlabel)
                                    % the state with highest index is the
                                    % artefact state
                                    nstatess = size(pe_by_state,4) - 1; % don't use the artefact state
                                else
                                    nstatess = size(pe_by_state,4);
                                end
                                z=spalloc(npes*nstatess,npes,npes*nstatess);
                                pe_tmp=zeros(npes*nstatess,1);
                                covpe_tmp=spalloc(npes*nstatess,npes*nstatess,npes*nstatess*nstatess);
                            end;

                            for iTime = 1 : size(pe_by_state,1), % indexes time points within trials

                                if(NKglm>1)

                                    % one state may be an artefact state
                                    if ~isempty(oat.source_recon.artefact_chanlabel)
                                        % the state with highest index is the
                                        % artefact state
                                        nstatess = size(pe_by_state,4) - 1; % don't use the artefact state
                                    else
                                        nstatess = size(pe_by_state,4);
                                    end

                                    for si = 1 : nstatess, % indexes state
                                        if(state_time{si}>=10),
                                            %z=[z; eye(npes)];
                                            z((si-1)*npes+1:(si)*npes,1:npes)=speye(npes);
                                            tmppe=permute(pe_by_state(iTime,:,f,si),[2 3 1]);

                                            if(strcmp(oat.first_level.cope_type,'coape')),
                                                tmppe=abs(tmppe);
                                            end;
                                            %pe_tmp=[pe_tmp; tmppe];
                                            pe_tmp((si-1)*npes+1:(si)*npes)=tmppe;

                                            tmpcovpe=permute(covarpe_by_state(iTime,:,:,f,si),[2 3 4 1]);
                                            %covpe_tmp=blkdiag(covpe_tmp,tmpcovpe);
                                            covpe_tmp((si-1)*npes+1:(si)*npes,(si-1)*npes+1:(si)*npes)=tmpcovpe;

                                            if(NKglm>1)
                                                first_level_results.cope_by_state(indind,iTime,:,f,si)=contrast_list_mat'*tmppe;
                                                first_level_results.stdcope_by_state(indind,iTime,:,f,si)=sqrt(diag(contrast_list_mat'*tmpcovpe*contrast_list_mat));
                                            end;
                                        end;
                                    end;

                                    [gam, covgam] = fixed_effects(pe_tmp,full(z),full(covpe_tmp));
                                    %gam=mean(cope_tmp);
                                else
                                    gam=permute(pe_by_state(iTime,:,f),[2 1]);
                                    covgam=permute(covarpe_by_state(iTime,:,:,f),[2 3 1]);
                                end;

                                if(NKglm>1)
                                    % if there is enough data,
                                    % average the artefact state over all
                                    % timepoints and conditions
                                    if ~isempty(oat.source_recon.artefact_chanlabel)
                                        if(state_time{nstatess+1}>=10)
                                            tmp = mean(pe_by_state(:,:,:,nstatess+1),2);
                                            artefact_map(indind,:) = squeeze(mean(tmp,1));
                                        else
                                        end
                                    end
                                end

                                if(strcmp(oat.first_level.cope_type,'coape')),
                                    gam=abs(gam);
                                end;

                                % NOW: apply contrasts
                                first_level_results.cope(indind,iTime,:,f)=contrast_list_mat'*gam;

                                if(strcmp(oat.first_level.cope_type,'acope')),
                                    first_level_results.cope(indind,iTime,:,f)=abs(first_level_results.cope(indind,iTime,:,f));
                                end;

                                first_level_results.stdcope(indind,iTime,:,f)=sqrt(diag(contrast_list_mat'*covgam*contrast_list_mat));

                                if(NKglm==1)
                                    first_level_results.cope_by_state(indind,iTime,:,f,1)=first_level_results.cope(indind,iTime,:,f);
                                    first_level_results.stdcope_by_state(indind,iTime,:,f,1)=first_level_results.stdcope(indind,iTime,:,f);
                                end;

                            end; % for iTime

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %% moving average within
                            %% the window specified
                            if ~isempty(first_level.post_glm_time_moving_av_win_size)

                                fsample=1/(D_times(2)-D_times(1));

                                windowSize = first_level.post_glm_time_moving_av_win_size*fsample;

                                for con=1:size(first_level_results.cope,3),
                                    tmp=first_level_results.cope(indind,:,con,f);

                                    first_level_results.cope(indind,:,con,f)=moving(tmp,round(windowSize));
                                    tmp=first_level_results.stdcope(indind,:,con,f).^2;

                                    first_level_results.stdcope(indind,:,con,f)=moving(tmp,round(windowSize)).^(0.5);

                                end;
                            end

                            %% do bc on the copes
                            if(sum(first_level.bc)>0),
                                baseline_time_indices=first_level_results.times<first_level.baseline_timespan(2) & first_level_results.times>first_level.baseline_timespan(1);

                                if isempty(find(baseline_time_indices))
                                    error(['There are no time points available in the specified baseline ' mat2str(first_level.baseline_timespan) ]);
                                end;

                                for c=1:length(contrast_list),
                                    if(first_level.bc(c))
                                        first_level_results.cope(indind,:,c,f)=first_level_results.cope(indind,:,c,f)-mean(first_level_results.cope(indind,baseline_time_indices,c,f),2);
                                    end;
                                end;

                            end;

                        if(0)
                            indind
                            iContrast=[2,4,5,6];iFreq=2;
                            figure(1);iFreq=2;plot(squeeze(first_level_results.cope(indind,:,iContrast,iFreq)))

                        end;
                    else
                        first_level_results.cope(indind,:,1,:)=dat(:,:);
                        first_level_results.stdcope(indind,:,1,:)=1;
                        if(NKglm>1)
                            first_level_results.cope_by_state(indind,:,1,:)=dat(:,:);
                            first_level_results.stdcope_by_state(indind,:,1,:)=1;
                        end;
                    end;
                else
                    warning('GLM not run: first_level.doGLM == 0');
                    results_fnames{sub} = ['GLM NOT RUN'];
                end;

            end % if first_level.doGLM
        end % for indind = 1:Nvoxels_out
        ft_progress('close');

    end; % for f=1:tf_num_freqs

    %%%%%%%%%%%%%%%%%%%
    %% store results
    first_level_results.tf_time_indices_into_D_times=tf_time_indices_into_D_times;

    first_level_results.x=first_level.x;
    first_level_results.contrasts=first_level.contrast;
    first_level_results.first_level_contrast_name=first_level.contrast_name;

    first_level_results.source_recon_results_fname=source_recon_results_fname;
    first_level_results.first_level=first_level;
    first_level_results.source_recon=source_recon_results.source_recon;
    first_level_results.level=1;
    %first_level_results.source_recon_results_fname=source_recon_results.fname;
    first_level_results.name=first_level.name;
    first_level_results.recon_method=source_recon_results.recon_method;
    try
        first_level_results.session_name=source_recon_results.session_name;
    catch,
        first_level_results.session_name=source_recon_results.subject_name;
    end;
    first_level_results.artefact_map = artefact_map;

    if first_level.doGLM
        first_level_results.fname=[ first_level_results.session_name '_' first_level_results.name ];
    end;

    %%%%%%%%%%%%%%%%%%
    %% combine grads if need to
    if strcmp(first_level_results.recon_method,'none') && ~strcmp(first_level.sensor_space_combine_planars,'dont_combine'),
        first_level_results=oat_stats_combine_grads(oat,first_level_results);
    end;

    %%%%%%%%%%%%%%%%%%%
    %% if in sensor space write out stats as SPM MEEG objects
    if first_level.doGLM
        if strcmp(source_recon_results.recon_method,'none'), % sensor space analysis

            S4=[];
            S4.oat=oat;
            S4.stats=first_level_results;
            [first_level_results.D_tstat, first_level_results.D_cope]=oat_save_spm_stats(S4);

        end;
    end;

    %%%%%%%%%%%%%%%%%%%
    %% save results to disk
    if first_level.save_trialwise_data
        if ~isempty(first_level.trialwise_directory)
            trlwise_fname=fullfile(first_level.trialwise_directory, [first_level_results.session_name '_trlwise']);
        else
            trlwise_fname=fullfile(oat.source_recon.dirname, [first_level_results.session_name '_trlwise']);
        end
        disp(['Saving trialwise data in file ' trlwise_fname]);
        save(trlwise_fname, 'trlwise', '-v7.3');
        clear trlwise
        first_level_results.trialwise_datafile=trlwise_fname;

    end % if first_level.save_data

    if first_level.doGLM,
        first_level_results.fname=[ first_level_results.session_name '_' first_level_results.name ];
        disp(['Saving statistics in file ' oat.source_recon.dirname '/' first_level_results.fname]);
        oat_save_results( oat, first_level_results );
        results_fnames{sub}=first_level_results.fname;
    end % if first_level.doGLM

    if (strcmp(source_recon_results.recon_method,'none') | isfield(first_level,'mni_coords')),
    else
        disp(['To create niftii files from this use a call to oat_save_nii_stats']);
    end;

    %%%%%%%%%%%%%%%%%%%
    %% generate source recon web report for this session
    if first_level.doGLM
        report = oat_first_level_stats_report(oat,first_level_results,'',first_level_results.report);
        first_level_results.report=osl_report_add_sub_report(first_level_results.report, report);

    end;
end;

%%%%%%%%%%%%%%%%%%%
%% generate first level web report
first_level_results.report=osl_report_write(first_level_results.report, oat.results.report);

end
