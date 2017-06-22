function [D,figs,S] = osl_africa(D,varargin)
    % AfRICA - ArteFact Rejection using Independent Component Analysis
    % performs ICA denoising of MEG data using either semi-manual or automated
    % identification artefact components.
    %
    % INPUTS:
    % - D           - SPM MEG object filename
    % - varargin    - key-value pairs, see inputParser below
    %
    % Romesh Abeysuriya 2017
    % Robert Becker 2017
    % Written by Henry Luckhoo and Adam Baker (pre-2017)

    arg = inputParser;
    arg.addParameter('modality','MEG'); % modality to use, default = 'MEG'
    arg.addParameter('montagename','AFRICA denoised data'); % New montage will be added with this name
    arg.addParameter('do_plots',false); % produce diagnostic plots, default = 0
    arg.addParameter('precompute_topos',true); % pre-compute and save IC spatial map topos after ica is computed for use in ident
    arg.addParameter('do_ica',~isfield(D,'ica')); % Do ICA decomposition step
    arg.addParameter('do_ident',true); % Do identification step
    arg.addParameter('do_remove',true); % Do removal step
    arg.addParameter('artefact_channels',{}); % Passed to ident_func
    arg.addParameter('ident_func',@identify_artefactual_components_manual);
    arg.addParameter('ident_params',struct); % Extra parameters for ident_func
    arg.addParameter('used_maxfilter',false); % Reduce ICA dimension if maxfilter was used, 
    arg.addParameter('ica_params',struct); % ICA parameters passed to run_sensorspace_ica
    arg.parse(varargin{:});
    S = arg.Results; % Result of parsing arguments is essentially the settings struct

    % Pass some of the options on to the identification function
    if ~isfield(S.ident_params,'artefact_channels')
        S.ident_params.artefact_channels = S.artefact_channels;
    end
    
    if ~isfield(S.ident_params,'modality')
        S.ident_params.modality = S.modality;
    end

    if ~isfield(S.ident_params,'do_plots')
        S.ident_params.do_plots = S.do_plots;
    end


    figs = struct('handles',[],'names',[],'titles',[]);

    D = D.montage('switch',0);

    % Do ICA decomposition
    if S.do_ica
        D = perform_sensorspace_ica(D,S); % Writes a new D.ica, also overwrites any old topos

        if S.precompute_topos
            topos = [];
            modalities = unique(D.chantype(find(strncmpi(S.modality,D.chantype,3)))); %#ok
            for m = 1:numel(modalities)
                disp(['Precomputing sensor topographies for modality ' modalities{m}]);
                topos = [topos component_topoplot(D,D.ica.sm,modalities(m))];
            end
            D.ica.topos   = topos;
        else
            D.ica.topos = []; % Compute in identify_artefactual_components_manual_gui.m
        end

        D.save(); % Consider taking this out - but ICA is so time consuming it might be worth it
    else
        assert(isfield(D,'ica'),'No precomputed results - D.ica is missing')
        fprintf('Using existing ICA decomposition\n')
    end

    % Identify bad components, store them together with metrics in the D object
    if S.do_ident
        [D.ica.bad_components, D.ica.metrics, figs] = S.ident_func(D,S.ident_params);
    else
        fprintf('Using existing bad_components\n')
    end

    % Remove bad components
    %%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS FROM THE DATA %%%%%%%%%%%%%%%%%%%%
    if S.do_remove
        D = remove_bad_components(D,S);
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

    % Good channels:
    chan_inds = indchantype(D,chantype,'GOOD');

    % Good timepoints/trials
    good_samples = ~all(badsamples(D,':',':',':'));
    good_samples = reshape(good_samples,1,D.nsamples*D.ntrials);

    % Select data:
    icadata = D(chan_inds,:,:);

    % Remove trial structure:
    icadata = reshape(icadata,size(icadata,1),[]);

    % Select good timepoints
    icadata = icadata(:,good_samples);

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

    if S.do_plots
        figure;
        semilogy(eigs_preNorm);
        ho;
        semilogy(eigs_postNorm,'r--');
        title('Raw and normalised eigen spectra'); legend('Raw', 'Normalised');
    end

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
    D.ica.chan_inds = chan_inds;
    D.ica.norm_vec = norm_vec;
    D.ica.sm = bsxfun(@times,sm,norm_vec);
    D.ica.topos = []; % Store topos later
    D.ica.metrics = []; % Store artefact metrics later
    D.ica.bad_components = []; % With new components, no bad components selected yet 

    % The commands below expand D.ica.sm out
    % Need to use D.ica.sm(chan_inds,:) elsewhere but that's what's expcted by a bunch of other things
    sm_full              = zeros(D.nchannels, size(sm,2));
    sm_full(chan_inds,:) = D.ica.sm; % So this is mapped onto all of the data w
    D.ica.sm = sm_full;

    % Reconstruct the IC timecourses with
    % tc = D(chan_inds,:,:)'*pinv(D.ica.sm)'
    %    or 
    % tc = D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))'
end


function D = remove_bad_components(D,S)
    % Take in a D object with D.ica.bad_components
    % Make a new montage with the components in D.ica.bad_components removed

    D = D.montage('switch',0);

    if strcmp(S.modality,'EEG')
        chantype = 'EEG';
        modality = 'EEG';
    else
        chantype = {'MEG','MEGANY'};
        modality = 'MEG';
    end

    % Good channels:
    chan_inds = indchantype(D,chantype, 'GOOD');
    tmp = struct(D);
    channels = tmp.channels(chan_inds);

    badchannels    = D.badchannels;
    bad_components = unique(D.ica.bad_components);
    megdata        = D(chan_inds,:,:);
    megdata        = reshape(megdata,size(megdata,1),[]);

    %%%%%%%%%%%%%%%%%%% REMOVE BAD COMPONENTS USING MONTAGE %%%%%%%%%%%%%%%%%%%

    sm = D.ica.sm;
    tc = (D(D.ica.chan_inds,:,:)'*pinv(D.ica.sm(D.ica.chan_inds,:))').';

    tra = eye(D.nchannels);
    dat_inv = pinv_plus(megdata', D.ica.params.num_ics);
    tra(chan_inds,chan_inds) = (eye(numel(chan_inds)) - dat_inv*(tc(bad_components,:)'*sm(chan_inds,bad_components)'))';
    labels = D.chanlabels;

    xchans = setdiff(1:D.nchannels,indchantype(D,chantype));
    tra(xchans,:)    = [];
    % tra(:,xchans)    = [];
    labels(xchans) = [];
   
    D = add_montage(D,tra,S.montagename,labels,channels)
   
end


