function [ results_fnames first_level_results ] = osl_run_first_level_continuous( oat )

% [ results_fnames ] = osl_run_first_level_continuous( oat )
%
% takes in an OAT, which needs to be setup by calling oat=osl_setup_oat(S), struct
% and runs first level continuous time GLM 
%
% This function should normally be called using osl_run_oat(oat);
%
% MWW 2012

OSLDIR = getenv('OSLDIR');

use_classes_in_glm=oat.first_level.hmm_do_glm_statewise;

first_level=oat.first_level;
source_recon_name='source_recon';

if strcmp(oat.source_recon.modalities{1},'EEG')
    modality_meeg='EEG';
else
    modality_meeg='MEG';
end

% check if user has asked not to do GLM
if(~isfield(first_level,'doGLM')),
    first_level.doGLM = 1;
end

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
%% check some of the settings

if(isfield(first_level,'connectivity_seed_mni_coord') && ~isempty(first_level.connectivity_seed_mni_coord)),
    error('This func does not handle seed connectivity analysis');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over sessions

clear stats;

% set first level diagnostic report up    
report_dir=[oat.results.plotsdir '/' oat.first_level.name];
first_level_results.report=osl_report_setup(report_dir,['First level (continuous)']);   

for subi_todo=1:length(first_level.sessions_to_do),   
        
    sub=first_level.sessions_to_do(subi_todo);
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(['%%%%%%%%%%%%%%%%%%%%%%%  RUNNING OAT FIRST LEVEL ON SESS = ' num2str(sub) '  %%%%%%%%%%%%%%%%%%%%%%%'])
        
    %sub/length(first_level.submatfiles),
    
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

            source_recon_results.is_sensor_space=1;
            
            chanindmeg = strmatch(modality_meeg, D.chantype);
            first_level_results.mask_indices_in_source_recon=1:length(chanindmeg);

            first_level_results.chanind=chanindmeg;
            first_level_results.chanlabels=D.chanlabels(chanindmeg);
            first_level_results.chantype=D.chantype(chanindmeg);
            
            first_level_results.mask_indices_in_source_recon=1:length(chanindmeg);
            
        else,

            source_recon_results.is_sensor_space=0;
            
            if(isfield(first_level,'mni_coords'))                         
                
                current_level_mni_coord = first_level.mni_coords;
                
                % load in lower level coords
                lower_level_mni_coord=source_recon_results.mni_coords;               
                %lower_level_mask_fname  = [oat.source_recon.dirname '/' source_recon_name '_mask'];
                %[ lower_level_mni_coord, ~ ] = osl_mnimask2mnicoords(lower_level_mask_fname);

                [~, iA, mask_indices]=intersect(current_level_mni_coord, lower_level_mni_coord,'rows');

                % sort indices to be in the order of the lower level mni_coords:
                [~, gg]=sort(iA);
                indices_in_lower_level=mask_indices(gg);

                first_level_results.mask_indices_in_source_recon = indices_in_lower_level;
                first_level_results.mni_coords = lower_level_mni_coord(indices_in_lower_level,:);
                
                % the mni coordinates from the specified list
                str=['Desired MNI coordinates ' ];
                for vox=1:size(current_level_mni_coord,1),
                    str=[str ', [' num2str(current_level_mni_coord(vox,:)) ']'];
                end;
                disp(str);
                
                % the mni coordinates from the lower level mask
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
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% set up potential time frequency bases in tf_settings, note that any time or
        % frequency averaging will be done AFTER the time-frequency decomposition
        
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
            if(freq_range(2)>post_tf_ds_factor*D.fsample/2)
                warning(['Nyquist=' num2str(post_tf_ds_factor*D.fsample/2) ', freq_range(2)>post_tf_ds_factor*D.fsample/2. Results will be dodgey due to looking at frequencies higher than the Nyquist limit.']);
            end;
        end;
        
        D_time_indices=find(source_recon_results.samples2use); % time indices in D.time used by source recon
        D_times = D.time(D_time_indices);
        
        if isempty(first_level.time_range)
            first_level.time_range=source_recon_results.woi;
        end
        
        tf_settings = [];
        tf_settings.tf_method               = first_level.tf_method;
        %tf_settings.tf_logtransform         = first_level.tf_logtransform;
        tf_settings.tf_logtransform         = 0;
        tf_settings.tf_hilbert_freq_res     = first_level.tf_hilbert_freq_res;
        tf_settings.tf_hilbert_freq_ranges  = first_level.tf_hilbert_freq_ranges;
        tf_settings.tf_freq_range           = freq_range; % set above
        tf_settings.tf_num_freqs            = first_level.tf_num_freqs;
        tf_settings.tf_multitaper_ncycles   = first_level.tf_multitaper_ncycles;
        tf_settings.raw_times               = D_times;
        tf_settings.ds_factor               = post_tf_ds_factor; 
        tf_settings.tf_morlet_factor        = first_level.tf_morlet_factor;
%        tf_settings.tf_multitaper_taper     = first_level.tf_multitaper_taper;
%        tf_settings.tf_multitaper_twin      = first_level.tf_multitaper_twin;
%        tf_settings.tf_multitaper_freqsmooth= first_level.tf_multitaper_freqsmooth;
        tf_settings.tf_calc_amplitude=0; % we will calc the amplitude after applying recon weights        
        tf_settings.tf_hilbert_do_bandpass_for_single_freq = first_level.tf_hilbert_do_bandpass_for_single_freq;
    
        out = osl_tf_transform( tf_settings ); % passed without a data matrix, osl_tf_transform just gets container sizes
        
        if isfield(out, 'tf_morlet_basis')
            first_level.tf_morlet_basis = out.tf_morlet_basis;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% setup time and frequency windows/indices for using later

        % setup result time dimension        
        tf_time_indices_into_D_times = find( ( out.tf_times >= first_level.time_range(1) ) & ( out.tf_times <= first_level.time_range(2) ));
        
%         tf_time_indices_into_D_times = false(1,numel(out.tf_times));
%         for i=1:size(first_level.time_range,1)
%              tf_time_indices_into_D_times(( out.tf_times >= first_level.time_range(i,1) ) & ( out.tf_times < first_level.time_range(i,2) )) = true;
%         end

        %tf_time_indices_into_D_times=find(tf_time_indices_into_D_times);
        tf_out_times=out.tf_times(tf_time_indices_into_D_times);        
        
        first_level_results.times=1;       
        if first_level.time_average
            error('Time averaging over whole window is not appropriate on continuous data');   
        end        
                
        % setup result freq dimension
        if ~isfield(out,'tf_freq_ranges')
            out.tf_freq_ranges=repmat(out.tf_freqs,2,1)';
        end;
        
        first_level_results.frequency_ranges=out.tf_freq_ranges;        
        first_level_results.frequencies = out.tf_freqs;
        
        % setup result spatial dimension
        Nvoxels_out=length(first_level_results.mask_indices_in_source_recon);

        clear out

        %%%%%%%%%%%%%%%%
        %% results containers
        ntpts                       = length(first_level_results.times);
        nfreqs                      = length(first_level_results.frequencies);
        
        if(first_level.doGLM)
            first_level_results.stdcope=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs]);% num_voxels x  num_timepoints x num_contrasts x num_freqs x Nvoxels_multidipole
            first_level_results.cope=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs]);% num_voxels x num_timepoints x num_contrasts x num_freqs x Nvoxels_multidipole
        end;

        first_level_results.pseudo_zstat_var = zeros([Nvoxels_out,1],'single');

    %end; % end of if(sub==1)

    disp(['Reconstruct time courses and computing stats for dataset ' source_recon_results.fname]);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% setup HMM 
    if strcmp(source_recon_results.recon_method,'none'),
        NKrecon=1;
        NKglm=NKrecon;
    else
        if strcmp(oat.source_recon.modalities{1},'EEG')
            NKrecon=numel(source_recon_results.BF.inverse.EEG.class);
        else
            NKrecon=numel(source_recon_results.BF.inverse.MEG.class);
        end

        NKglm=NKrecon;

        if(~use_classes_in_glm),
            NKglm=1;
        end;

        if(NKglm>1)
            if( post_movingaverage_ds_factor~=1 ),
                error('Can not do any downsampling when using multiple classes (e.g. HMM)');
            end;
        end;
    end;
    
    if(first_level.doGLM)
        
        first_level_results.stdcope_by_state=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs,NKglm]);% num_voxels x  num_timepoints x num_contrasts x num_freqs x Nvoxels_multidipole
        first_level_results.cope_by_state=zeros([Nvoxels_out,ntpts,length(first_level.contrast),nfreqs,NKglm]);% num_voxels x num_timepoints x num_contrasts x num_freqs x Nvoxels_multidipole
    end;
       
    %%%%%%%%%%%%%%%%
    %% setup the GLM design matrix
    % design matrix should be num_PEs x num_timepoints           
    if(~isfield(first_level,'design_matrix')),

        if(first_level.doGLM)
            error('Design matrix needs to be specified');
        else
            % setup dummy design matrix
            first_level.design_matrix=ones(1,length(D_time_indices));
        end;
    end;

    if subi_todo == 1
        if  strcmp(first_level.design_matrix,'avg')
            avg_design_matrix = true;
        else
            avg_design_matrix = false;
        end
    end

    if avg_design_matrix == true
        disp('Using single constant regressor as design matrix');
        
        if post_tf_ds_factor~=1
            first_level.design_matrix = ones(1,numel(tf_out_times));       
        else
            first_level.design_matrix = ones(1,length(D.time));
            % first_level.design_matrix should have same ntpts as D_time
            if(size(first_level.design_matrix,2)~=length(D.time)),
                error(['Design matrix dimensions [npes x ntpts] = ' num2str(size(first_level.design_matrix)) ' do not match the data: ntpts=' num2str(length(D_time_indices))]); 
            end;
        end
        processed_design_matrix = first_level.design_matrix;
    else
        
        % first_level.design_matrix should have same ntpts as D_time
        if(size(first_level.design_matrix,2)~=length(D.time)),
            error(['Design matrix dimensions [npes x ntpts] = ' num2str(size(first_level.design_matrix)) ' do not match the data: ntpts=' num2str(length(D_time_indices))]); 
        end;
        processed_design_matrix=first_level.design_matrix(:,D_time_indices);    
        
        % temporally downsample to mimic what happens to data 
        if post_tf_ds_factor~=1,
            processed_design_matrix = osl_spm_resample(processed_design_matrix,post_tf_ds_factor);  
        end;
    end

    disp(['Design matrix dimensions [npes x ntps] = ' num2str(size(first_level.design_matrix))]);

    if(length(oat.first_level.bc)~=length(oat.first_level.contrast)),
        error('first_level.bc and first_level.contrasts need to be the same length');
    end;
    
    if(length(oat.first_level.contrast_name)~=length(oat.first_level.contrast))
        error('oat.first_level.contrast and oat.first_level.contrast_name need to be same length');
    end;
    
    % apply moving average to design matrix
    if ~isempty(first_level.time_moving_av_win_size)
        fsample=1/(D_times(2)-D_times(1));        
        windowSize = first_level.time_moving_av_win_size*fsample;
                                
        for pp=1:size(processed_design_matrix,1),
            processed_design_matrix(pp,:)=moving(processed_design_matrix(pp,:),round(windowSize));
        end;
    end;

    processed_design_matrix2 = zeros(size(processed_design_matrix,1),length(tf_time_indices_into_D_times));
    for pp=1:size(processed_design_matrix,1),
        processed_design_matrix2(pp,:)=processed_design_matrix(pp,tf_time_indices_into_D_times);
    end;   
    processed_design_matrix=processed_design_matrix2;
    
    disp(['Design matrix dimensions [npes x ntpts] = ' num2str(size(processed_design_matrix))]);
        
    % contrasts                                  
    contrast_list=first_level.contrast;
    for c=1:length(contrast_list),
        contrast_list{c}=contrast_list{c}(:);
    end;
    
    Ntrials=1;
        
    triallist=1;
        
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
    if(~use_classes_in_glm),
        class_samples_inds_glm{1}=ones(1,length(D_time_indices),triallist);
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
    %% Do T-F transform in sensor space
    %% Transform the data - either downsampling in time domain, or a TF transform
    
    disp('Doing T-F transform in sensor space...');

    tf_settings = [];
    tf_settings.tf_method               = first_level.tf_method;
    %tf_settings.tf_logtransform         = first_level.tf_logtransform;
    tf_settings.tf_logtransform         = 0;    
    tf_settings.tf_hilbert_freq_res     = first_level.tf_hilbert_freq_res;
    tf_settings.tf_freq_range           = freq_range; % set above
    tf_settings.tf_hilbert_freq_ranges  = first_level.tf_hilbert_freq_ranges;
    tf_settings.tf_num_freqs            = first_level.tf_num_freqs;
    tf_settings.raw_times               = D_times;
    tf_settings.ds_factor               = post_tf_ds_factor; 
    tf_settings.tf_morlet_factor        = first_level.tf_morlet_factor;
%    tf_settings.tf_multitaper_taper                = first_level.tf_multitaper_taper;
%    tf_settings.tf_multitaper_twin                 = first_level.tf_multitaper_twin;
%    tf_settings.tf_multitaper_freqsmooth           = first_level.tf_multitaper_freqsmooth;
    tf_settings.tf_calc_amplitude=0; % we will calc the amplitude after applying recon weights                
    tf_settings.tf_hilbert_do_bandpass_for_single_freq = first_level.tf_hilbert_do_bandpass_for_single_freq;
    
    if isfield(first_level,'tf_morlet_basis')
        tf_settings.tf_morlet_basis = first_level.tf_morlet_basis;
    end
    
    sensor_data=D(:, D_time_indices, triallist);
    
    dims=size(sensor_data);
    sensor_data_tf=zeros(dims(1),length(tf_time_indices_into_D_times),Ntrials,first_level.tf_num_freqs);
        
    if(first_level.tf_num_freqs>1), ft_progress('init', 'etf');end;    
    for indind=1:size(sensor_data,1),
        if(first_level.tf_num_freqs>1), ft_progress(indind/size(sensor_data,1));end;

        for tri=1:size(sensor_data,3),
            out = osl_tf_transform( tf_settings , sensor_data(indind,:,tri) ); % do the transformation
            %sensor_data_tf(indind,:,tri,:)=permute(out.dattf,[2 3 1]);
            sensor_data_tf(indind,:,tri,:)=out.dattf;
        end;
    end;
    class_samples_inds_glm_tf=class_samples_inds_glm;
    class_samples_inds_recon_tf=class_samples_inds_recon;
    
    if(first_level.tf_num_freqs>1), ft_progress('close');end;
      
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
    %
    % We need to be careful to preseve both the real and imaginary parts of
    % complex data if necessary. They are stored in separate files for the
    % moment and will be separately reconstructed from the beamformer
    % weights.

    Sc=[];
    Sc.D = D;
    [path nm ext]=fileparts(fullfile(D));
    Sc.newname = [path '/TF' nm ext];
    Sc.newdata = real(sensor_data_tf);
    Sc.time = tf_out_times;
    Sc.frequencies = first_level_results.frequencies;
    Sc.remove_montages=0;
    D_tf = osl_change_spm_eeg_data( Sc );   

    % add back in Class channel:        
    classchanind=find(strcmp(D.chanlabels,'Class'));
    classchanind_tf=find(strcmp(D_tf.chanlabels,'Class'));
    if D_tf.nfrequencies==1
        D_tf(classchanind_tf,:,1)=D(classchanind,D_time_indices(tf_time_indices_into_D_times),1);
    else
        for ff = 1:D_tf.nfrequencies
            D_tf(classchanind_tf,ff,:,1)=permute(D(classchanind,D_time_indices(tf_time_indices_into_D_times),1),[1 3 2]);
        end
    end
    
    if ~isreal(sensor_data_tf)
        Sc.newname = [path '/TFimag' nm ext];
        Sc.newdata = imag(sensor_data_tf);
        D_tf_imag = osl_change_spm_eeg_data( Sc );

        if D_tf_imag.nfrequencies==1
            D_tf_imag(classchanind_tf,:,1)=D(classchanind,D_time_indices(tf_time_indices_into_D_times),1);
        else
            for ff = 1:D_tf_imag.nfrequencies
                D_tf_imag(classchanind_tf,ff,:,1)=permute(D(classchanind,D_time_indices(tf_time_indices_into_D_times),1),[1 3 2]);
            end
        end
    else
        D_tf_imag = [];
    end
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

    %%%%%%%%%%%%%%%%
    %% Setup design matrices that can be reused over voxels
    
    const_reg=0;        
    for p=1:size(processed_design_matrix,1),
        if(range(processed_design_matrix(p,:))==0),
            if(const_reg),
                warning('Multiple constant regressors in the design matrix!!!!?');
            end;
            const_reg=1;

        end;                                    
    end;

    if(const_reg && oat.first_level.do_glm_demean)
        warning('CAREFUL: are you sure you want the first_level.do_glm_demean flag on with constant regressors in the design matrix!!!!?');
    elseif(~const_reg && ~oat.first_level.do_glm_demean)
        warning('CAREFUL: are you sure you want no first_level.do_glm_demean flag on with no constant regressors in the design matrix!!!!?');
    end;       
      
    ntpts=1;
    x_stored=cell(NKglm,ntpts);
    pinvxtx_stored=cell(NKglm,ntpts);
    pinvx_stored=cell(NKglm,ntpts);
    
    for si=1:NKglm,
        class_samples_inds_tf{si}=class_samples_inds_glm_tf{si}(1,:,:);
        state_tinds=class_samples_inds_tf{si};
       
        total_time=(D_times(end)-D_times(1))*size(state_tinds,3);
        state_time=total_time*sum(sum(squeeze(state_tinds)))/prod(size(squeeze(state_tinds)));
 
        if(state_time>10),
            x=processed_design_matrix';
            
            x=x(find(permute(state_tinds(1,:),[2 1])),:);                        

            if(oat.first_level.do_glm_demean)
                for pp=1:size(x,2),
                    x(:,pp)=x(:,pp)-mean(x(:,pp));
                end;
            end;

            x_stored{si}=x;
            pinvxtx_stored{si}=pinv(x'*x);
            pinvx_stored{si}=pinv(x);
        end;
    end;
    clear state_tinds;
    
    contrast_list_mat=cell2mat(contrast_list);
    
    %%%%%%%%%%%%%%%%
    %% do everything from now on one voxel at a time (for the sake of
    %% RAM) in the first level mask space                
    if(first_level.doGLM)
        disp('First level COPEs outputted will have dimension Nvoxels x Ntpts x Ncontrasts x Nfreqs:');
        tmp=[size(first_level_results.cope) 1 1];
        disp([num2str(tmp(1:4))]);
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now loop over voxels
    ft_progress('init', 'etf');
    
    % debug:
    %[indfind, mni_coord]=nearest_vec(source_recon_results.mni_coords,[-32 -26 50]);
      
    S2=[];
    S2.class_samples_inds=class_samples_inds_recon_tf;
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
    
    for indind=1:Nvoxels_out, % indexes brain space      
        
        ft_progress(indind/Nvoxels_out);                                   

        S2.index=first_level_results.mask_indices_in_source_recon(indind);     

        if source_recon_results.is_sensor_space==1
            dat_tf = permute(S2.D(indind,:,:,1),[3 2 1 4]); %[1 3 4 2]
        else
            [dat_tf S2] = osl_get_recon_timecourse( S2 );
            if isa(D_tf_imag,'meeg')
                Sc.D = D_tf_imag;
                [tmp ~] = osl_get_recon_timecourse( S2 );
                dat_tf = dat_tf + 1i.*tmp;
                clear tmp
            end
            % dat_tf needs to trials x time:
            if isempty(S2.D.nfrequencies)
                dat_tf=permute(dat_tf,[2 1 3]); %[2 1 3]
            else
                dat_tf=permute(dat_tf,[3 2 1]);
            end
        end
        
        first_level_results.pseudo_zstat_var(indind)=var(dat_tf(:,1));

        if sum(isnan(squash(dat_tf)))~=0,

            if(first_level.doGLM)
                first_level_results.cope(indind,:,:,:)=nan;
                first_level_results.stdcope(indind,:,:,:)=nan;         
                first_level_results.cope_by_state(indind,:,:,:,:)=nan;
                first_level_results.stdcope_by_state(indind,:,:,:,:)=nan; 
            end;
            error('NANs in recon data');

        else
               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %% calc the amplitude (up to now dat_tf is complex)
            dat_tf=abs(dat_tf);
            
            dat=dat_tf; % num_tpts x num_freqs
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %% moving average within
            %% the window specified
            if ~isempty(first_level.time_moving_av_win_size)

                fsample=1/(D_times(2)-D_times(1));        
                windowSize = first_level.time_moving_av_win_size*fsample;
                
                clear dat_new;
                
                for f=1:size(dat,2),
                    tmp=permute(dat(:,f),[1,2,3]);
                    tmp=moving(tmp,round(windowSize));
                    dat_new(:,f)=tmp;
                end;

                dat=dat_new;                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %% temporally downsample (needs to be done after abs (for t-f data))                                
            if post_movingaverage_ds_factor~=1,
                clear dat_new;
                for f=1:size(dat,2),
                    tmp=permute(dat(:,f),[1,2,3]);
                    tmp=osl_spm_resample(tmp',post_movingaverage_ds_factor)';  
                    dat_new(:,f)=tmp;
                end;
                dat=dat_new;
            end;
            
            glm_input_times=linspace(first_level.time_range(1),first_level.time_range(2),length(dat));
                                                                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            %% now do the GLM stats
            if first_level.doGLM
                npes=size(contrast_list_mat,1);
                
                if(size(dat,1)>1)                       
                    for f = 1 : size(dat,2),
                        clear total_time state_time;
                        for si=1:NKglm,

                            state_tinds=class_samples_inds_glm_tf{si};

                            total_time=(D_times(end)-D_times(1))*size(state_tinds,3);
                            state_time{si}=total_time*sum(sum(squeeze(state_tinds)))/prod(size(squeeze(state_tinds)));

                            if(state_time{si}<10),

                                if(NKglm>1)
                                    %first_level_results.cope_by_state(indind,:,:,f,si)=0;
                                    %first_level_results.stdcope_by_state(indind,:,:,f,si)=inf;

                                    pe_by_state(:,:,si)=0;
                                    for t= 1 : size(dat,2),
                                        covarpe_by_state(t,:,:,si)=inf;
                                    end;

                                end;
                            else

                                t=1;
                                x=x_stored{si};
                                pinvxtx=pinvxtx_stored{si};
                                pinvx=pinvx_stored{si};                               

                                dattf=dat(find(permute(state_tinds(1,:),[2 1])),f,t);

                                if(oat.first_level.do_glm_demean)
                                    dattf=demean(dattf);                                                        
                                end;
                                                                        
                                [peout, covarpeout]=glm_fast_for_meg(dattf,x,pinvxtx,pinvx);

                                pe_by_state(:,:,si)=reshape(peout,npes,size(dat,3))'; 

                                covarpe_by_state(t,:,:,si)=covarpeout((t-1)*npes+1:t*npes,(t-1)*npes+1:t*npes);

                            end;                    
                                                   
                        end; % si (state index)

                        %% do fixed effects averaging of states
                        for iTime = 1 : size(pe_by_state,1), % indexes time points within trials

                            if(NKglm>1)
                                pe_tmp=[];
                                covpe_tmp=[];
                                z=[];

                                for si = 1 : size(pe_by_state,3), % indexes state
                                    if(state_time{si}>=10),
                                        z=[z; eye(npes)];
                                        tmppe=permute(pe_by_state(iTime,:,si),[2 3 1]);
                                        pe_tmp=[pe_tmp; tmppe];
                                        
                                        if(strcmp(oat.first_level.cope_type,'coape')),
                                            tmppe=abs(tmppe);
                                        end;
                                        first_level_results.cope_by_state(indind,iTime,:,f,si)=contrast_list_mat'*tmppe;

                                        tmpcovpe=permute(covarpe_by_state(iTime,:,:,si),[2 3 4 1]);
                                        covpe_tmp=blkdiag(covpe_tmp,tmpcovpe);
                                        
                                        first_level_results.stdcope_by_state(indind,iTime,:,f,si)=sqrt(diag(contrast_list_mat'*tmpcovpe*contrast_list_mat));                                                                
                                    end;                                                                
                                end;

                                [gam, covgam] = fixed_effects(pe_tmp,z,covpe_tmp);                                
                                
                                %gam=mean(cope_tmp);
                            else
                                gam=permute(pe_by_state(iTime,:),[2 1]);
                                covgam=permute(covarpe_by_state(iTime,:,:),[2 3 1]);                                
                            end;

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
                                               
                        %% do bc on the copes
                        if(sum(first_level.bc)>0),
                            warning('No baseline correction is carried out with continuous time time-wise GLM');                            
                        end;
                    end;
                    
                else
                    first_level_results.cope(indind,:,1,:)=dat(:,:);
                    first_level_results.stdcope(indind,:,1,:)=1;                    
                    first_level_results.cope_by_state(indind,:,1,:)=dat(:,:);
                    first_level_results.stdcope_by_state(indind,:,1,:)=1;
                end;
            
            else,
                if(~isfield(first_level_results,'glm_input_data')),
                    first_level_results.glm_input_times=glm_input_times;
                    first_level_results.glm_input_frequencies=out.tf_freqs;
                    first_level_results.glm_input_data=zeros([Nvoxels_out, size(dat)]);% num_voxels x  num_timepoints x num_contrasts x num_freqs                    
                end;                    
                    
                first_level_results.glm_input_data(indind,:,:,:)=dat;
            end; % if first_level.doGLM                   
        end;
    end; % for indind = 1:Nvoxels_out
              
    ft_progress('close');    
            
    
    %%%%%%%%%%%%%%%%%%%
    %% store results
    
    first_level_results.tf_time_indices_into_D_times=tf_time_indices_into_D_times;
    first_level_results.source_recon_results_fname=source_recon_results_fname;
    first_level_results.first_level=first_level;
    first_level_results.source_recon=source_recon_results.source_recon;
    first_level_results.first_level_contrast_name=first_level.contrast_name;
    
    first_level_results.level=1;
    %first_level_results.source_recon_results_fname=source_recon_results.fname;
    first_level_results.name=first_level.name;
    first_level_results.recon_method=source_recon_results.recon_method;
    first_level_results.session_name=source_recon_results.session_name;
          
    first_level_results.fname=[ first_level_results.session_name '_' first_level_results.name ];
            
    if first_level.doGLM
        first_level_results.x=x;
        first_level_results.contrasts=first_level.contrast;
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
    %% save results
    
    disp(['Saving statistics in file ' oat.source_recon.dirname '/' first_level_results.fname]);
    oat_save_results( oat, first_level_results );
    results_fnames{sub}=first_level_results.fname;

    if (strcmp(source_recon_results.recon_method,'none')),
    else
        disp(['To create niftii files from this use a call to oat_save_nii_stats']);
    end;

    %%%%%%%%%%%%%%%%%%%
    %% generate source recon web report for this session
    if first_level.doGLM
        report = oat_first_level_stats_report(oat,first_level_results,'',first_level_results.report);
        first_level_results.report=osl_report_add_sub_report(first_level_results.report, report)
    end;

end;

%%%%%%%%%%%%%%%%%%%
%% generate first level web report
first_level_results.report=osl_report_write(first_level_results.report, oat.results.report);        

end
