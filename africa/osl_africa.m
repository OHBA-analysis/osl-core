function D = osl_africa(D,varargin)
    % AfRICA - ArteFact Rejection using Independent Component Analysis
    % performs ICA denoising of MEG data using either semi-manual or automated
    % identification artefact components.
    %
    % INPUTS:
    % - D           - SPM MEG object
    % - varargin    - key-value pairs, see inputParser below
    %
    % Romesh Abeysuriya 2017
    % Robert Becker 2017
    % Written by Henry Luckhoo and Adam Baker (pre-2017)

    assert(strcmp(D.type,'continuous'),'osl_africa is only compatible with continuous MEEG objects'); % See issue #43 on GitHub

    arg = inputParser;

    % INPUT SETTINGS
    arg.addParameter('modality','MEG'); % modality to use, default = 'MEG'
    arg.addParameter('artefact_channels',{'EOG','ECG'},@iscell); % Specify which channels are artefactual e.g. {'EOG1','ECG',32}. Can be chantype OR channel index number
    arg.addParameter('mains_frequency',50); 
    arg.addParameter('used_maxfilter',false); % Reduce ICA dimension if maxfilter was used, 

    % STAGE SETTINGS
    arg.addParameter('do_ica',~isfield(D,'ica')); % Do ICA decomposition step
    arg.addParameter('do_ident','auto'); % Do identification step - options are: false/empty to skip,'auto' or 'manual'
    arg.addParameter('do_remove',true); % Do removal step

    %  ICA SETTINGS
    arg.addParameter('precompute_topos',true); % pre-compute and save IC spatial map topos after ica is computed for use in ident
    arg.addParameter('ica_params',struct,@isstruct); % ICA parameters passed to run_sensorspace_ica - typically do not require changing

    % AUTOMATIC IDENT SETTINGS
    arg.addParameter('auto_max_num_artefact_comps',10); % Maximum number of new components to reject for each reason
    arg.addParameter('auto_do_mains',false); % Used by manual and auto
    arg.addParameter('auto_mains_kurt_thresh',0.4);
    arg.addParameter('auto_do_kurt',false); % Used by manual and auto
    arg.addParameter('auto_kurtosis_thresh',20); 
    arg.addParameter('auto_kurtosis_wthresh',0); 
    arg.addParameter('auto_artefact_chans_corr_thresh',0.5,@isscalar);

    % OUTPUT SETTINGS
    arg.addParameter('montagename','AFRICA denoised data'); % New montage will be added with this name

    arg.parse(varargin{:});
    S = arg.Results; % Result of parsing arguments is essentially the settings struct

    original_montage_index = D.montage('getindex'); % If do_remove=False then the returned MEEG will have the original montage
    D = D.montage('switch',0);

    % Do ICA decomposition
    if S.do_ica
        D = perform_sensorspace_ica(D,S); % Writes a new D.ica, also overwrites any old topos

        if S.precompute_topos
            topos = [];
            for m = 1:numel(D.ica.modalities)
                disp(['Precomputing sensor topographies for modality ' D.ica.modalities{m}]);
                topos = [topos component_topoplot(D,D.ica.sm,D.ica.modalities(m))];
            end
            D.ica.topos   = topos;
        else
            D.ica.topos = [];
        end

        % In general, we avoid saving the D object so that users are
        % explicitly aware when changes are made on-disk. However, ICA is very
        % time consuming, so on balance this is more user-friendly
        fprintf(1,'** Saving changes to disk **\n');
        D.save(); 
    else
        assert(isfield(D,'ica'),'Skip running ICA was specified, but there are no precomputed results - D.ica is missing')
        fprintf('Using existing ICA decomposition\n')
    end

    % Identify bad components, store them together with metrics in the D object
    if S.do_ident
        [D.ica.metrics,tc] = compute_metrics(D,S.mains_frequency,S.artefact_channels);
        switch S.do_ident
            case 'auto'
                [D.ica.bad_components, D.ica.auto_reason] = identify_artefactual_components_auto(D,S);
            case 'manual'
                D.ica.bad_components = identify_artefactual_components_manual(D,tc,D.ica.topos,D.ica.metrics,D.ica.bad_components);
                for j = 1:length(D.ica.bad_components)
                    fprintf('IC %d marked bad\n',D.ica.bad_components(j));
                end
            otherwise
                error('Did not recognize ident type - must be auto, manual, or false/empty');
        end
    else
        fprintf('Using existing bad_components\n')
    end

    % Remove bad components
    if S.do_remove
        D = remove_bad_components(D,S);
    else
        D = D.montage('switch',original_montage_index);
    end

end % MAIN FUNCTION

function D = perform_sensorspace_ica(D,S)
    % Given a D object and a settings struct, run fastica
    % Return a D object with D.ica created or overwritten

    D = D.montage('switch',0);
        
    if strcmp(S.modality,'EEG')
        chantype = 'EEG';
    else
        chantype = {'MEG','MEGANY'}; % need both flags to capture MEG, MEGMAG, MEGGRAD and MEGPLANAR channel types.
    end

    % Good channels and timepoints/trials
    chan_inds = indchantype(D,chantype,'GOOD');
    ica_good_samples = good_samples(D,chan_inds);
    ica_good_samples = reshape(ica_good_samples,1,D.nsamples*D.ntrials);

    % Select data:
    icadata = D(chan_inds,:,:);

    % Remove trial structure:
    icadata = reshape(icadata,size(icadata,1),[]);

    % Select good timepoints
    icadata = icadata(:,ica_good_samples);

    %%%%%%%%%%%%%%%%%%%% APPLY MAXFILTER SPECIFIC SETTINGS %%%%%%%%%%%%%%%%%%%%
    if isfield(S,'used_maxfilter') && S.used_maxfilter
        num_ics_default = 62;
        mag_cutoff      = 62;
        plan_cutoff     = 62;
    else
        num_ics_default = 150;
        mag_cutoff  = sum(strcmp(D.chantype(chan_inds),'MEGMAG'))    - 5;
        plan_cutoff = sum(strcmp(D.chantype(chan_inds),'MEGPLANAR')) - 5;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% SET FASTICA PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%
    ica_params = inputParser;
    ica_params.addParameter('num_ics',num_ics_default);       
    ica_params.addParameter('last_eig',num_ics_default);      
    ica_params.addParameter('nonlinearity','tanh');  
    ica_params.addParameter('approach','symm');      
    ica_params.addParameter('stabilization','on'); 
    ica_params.addParameter('max_iter',1000);
    ica_params.parse(S.ica_params);
    ica_params = ica_params.Results;
    
    ica_params.num_ics  = min(ica_params.num_ics, size(icadata,1)); % added by DM
    ica_params.last_eig = min(ica_params.last_eig,size(icadata,1)); % added by DM

    %%%%%%%%%%%%%%%%%%%%  MINIMUM EIGENVALUE NORMALISATION %%%%%%%%%%%%%%%%%%%%
    if strcmp(S.modality,'EEG')  % added by DM
        norm_vec = max(abs(icadata(:)))/1000*ones(size(icadata,1),1);
    else
        norm_vec = ones(numel(chan_inds),1);
        if any(strcmp(D.chantype,'MEGMAG')) && any(strcmp(D.chantype,'MEGPLANAR'))
            mag_min_eig = svd(cov(icadata(strcmp(D.chantype(chan_inds),'MEGMAG'),:)'));
            mag_min_eig = mean(mag_min_eig(mag_cutoff-2:mag_cutoff));
            
            plan_min_eig = svd(cov(icadata(strcmp(D.chantype(chan_inds),'MEGPLANAR'),:)'));
            plan_min_eig = mean(plan_min_eig(plan_cutoff-2:plan_cutoff));
            
            norm_vec(strcmp(D.chantype(chan_inds),'MEGMAG'))    = mag_min_eig;
            norm_vec(strcmp(D.chantype(chan_inds),'MEGPLANAR')) = plan_min_eig;
        else
            norm_vec = norm_vec*min(svd(cov(icadata(:,:)')));
        end
        norm_vec = sqrt(norm_vec);
        
    end

    eigs_preNorm  = svd(cov(icadata'));

    % Apply normalisation
    icadata = icadata ./ repmat(norm_vec,1,size(icadata,2));

    eigs_postNorm = svd(cov(icadata'));


    %%%%%%%%%%%%%%%%%% AUTOMATIC DIMENSIONALITY ESTIMATION %%%%%%%%%%%%%%%%%%%%
    if 0 == ica_params.num_ics
        ica_params.num_ics = spm_pca_order(icadata');
    end
    if 0 == ica_params.last_eig
        ica_params.last_eig = ica_params.num_ics;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% ICA DECOMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [tc,sm,~] = fastica(icadata,                           ...
        'g',                ica_params.nonlinearity,  ...
        'lastEig',          ica_params.last_eig,      ...
        'numOfIC',          ica_params.num_ics,       ...
        'approach',         ica_params.approach,  ...
        'stabilization',    ica_params.stabilization, ...
        'maxNumIterations', ica_params.max_iter); % changed by DM

    if ica_params.num_ics ~= size(tc,1)
        fprintf('\n%s%d%s%d%s\n%','Data dimensionality insufficient to support ', ica_params.num_ics, ' components. Number of components has been reduced to ', size(tc,1), '.');
        ica_params.num_ics = size(tc,1);
    end

    D.ica = struct;
    D.ica.params = ica_params;
    D.ica.eigs_preNorm = eigs_preNorm;
    D.ica.eigs_postNorm = eigs_postNorm;
    D.ica.good_samples = ica_good_samples; % These were the good timepoints that went into the ICA. i.e. use these when indexing the ICs
    D.ica.chan_inds = chan_inds;
    D.ica.norm_vec = norm_vec;
    D.ica.sm = bsxfun(@times,sm,norm_vec);
    D.ica.topos = []; % Store topos later
    D.ica.metrics = []; % Store artefact metrics later
    D.ica.bad_components = []; % With new components, no bad components selected yet 
    D.ica.auto_reason = {}; % Record reason for automatic rejection (if used)
    D.ica.modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3)))); 

    % The commands below expand D.ica.sm out to a full montage
    sm_full              = zeros(D.nchannels, size(sm,2));
    sm_full(chan_inds,:) = D.ica.sm; % So this is mapped onto all of the data w
    D.ica.sm = sm_full;

    % Reconstruct the IC timecourses with
    % tc = D(:,:,:)'*pinv(D.ica.sm)'
    %    or 
    % tc = D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))'
end


function D = remove_bad_components(D,S)
    % Take in a D object with D.ica.bad_components
    % Make a new montage where all the channels stay the same but 
    % they have had the bad ICA components subtracted from them

    D = D.montage('switch',0);

    if strcmp(S.modality,'EEG')
        chantype = 'EEG';
        modality = 'EEG';
    else
        chantype = {'MEG','MEGANY'};
        modality = 'MEG';
    end

    chan_inds = D.ica.chan_inds; % This is a record of which channels have had components subtracted from them
    bad_components = unique(D.ica.bad_components);
    megdata        = D(chan_inds,:,:);
    megdata        = reshape(megdata,size(megdata,1),[]);

    tc = (D(chan_inds,:,:)'*pinv(D.ica.sm(chan_inds,:))').';
    tra = eye(D.nchannels);
    dat_inv = pinv_plus(megdata', D.ica.params.num_ics);
    tra(chan_inds,chan_inds) = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*D.ica.sm(chan_inds,bad_components)'))';
    
    tmp = struct(D);
    D = add_montage(D,tra,S.montagename,D.chanlabels,tmp.channels);
end


