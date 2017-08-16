% osl_ica_maps_group_stats
%
% Syntax: oil = osl_ica_maps_group_stats(oil);
%
% This is currently in developement.
%
%
%
% Henry Luckhoo (henry.luckhoo@trinity.ox.ac.uk)
% Version 1.3
% 25.03.13

function oil = oil_ica_maps_group_stats(oil)

OSLDIR = getenv('OSLDIR');
fprintf('Oil group-level analysis. \n');

%% Set-Up Variables
if isfield(oil.ica_group_level, 'group_design_matrix');                design_matrix                     = oil.ica_group_level.group_design_matrix;                             else, error('Group Design Matrix unspecified');      end
if isfield(oil.ica_group_level, 'group_contrast');                     contrast                          = oil.ica_group_level.group_contrast;                                  else, error('Contrasts missing');                    end
if isfield(oil.ica_group_level, 'comps2use');                          comps2use                         = oil.ica_group_level.comps2use;                                       else, comps2use = 1:oil.ica.num_ics;                 end
if isfield(oil.ica_group_level, 'Npermutations');                      nperms                            = oil.ica_group_level.Npermutations;                                   else, nperms    = 5000;                              end
if isfield(oil.ica_group_level, 'use_randomise'),                      use_randomise                     = oil.ica_group_level.use_randomise;                                   else, use_randomise = 1;                             end
if isfield(oil.ica_group_level, 'demean_designmatrix'),                demean_designmatrix               = oil.ica_group_level.demean_designmatrix;                             else, demean_designmatrix = 0;                       end
if isfield(oil.ica_group_level, 'group_varcope_spatial_smooth_fwhm'),  group_varcope_spatial_smooth_fwhm = oil.ica_group_level.group_varcope_spatial_smooth_fwhm;               else, group_varcope_spatial_smooth_fwhm = 0;         end
if isfield(oil.ica_group_level, 'cluster_threshold'),                  thresh                            = oil.ica_group_level.cluster_threshold;                               else, thresh = 2.3;                                  end
if isfield(oil.ica_first_level.results, 'cope_files');                 cope_files                        = oil.ica_first_level.results.cope_files;                              else, error('Single subject cope maps are missing'); end
if isfield(oil.ica_first_level.results, 'tstat_files');                tstat_files                       = oil.ica_first_level.results.tstat_files;                             else, error('Single subject t-stats are missing');   end
if isfield(oil.enveloping, 'gridstep'),                                gridstep                          = oil.enveloping.gridstep;                                             else, error('Gridstep Missing');                     end

mask_fname = [OSLDIR '/std_masks/MNI152_T1_' num2str(gridstep) 'mm_brain.nii.gz'];    
Nsubs       = length(cope_files);
Nvoxels     = size(nii.quickread(cope_files{1},gridstep),1);
Nics        = numel(comps2use);
Ncontrasts = length(contrast);

stats_dir = [oil.source_recon.dirname '/' oil.enveloping.name '/' oil.concat_subs.name '/' oil.ica.name '/' oil.ica_first_level.name '/' oil.ica_group_level.name];
if ~isdir(stats_dir), mkdir(stats_dir); end

%% Load in single subject copes
copes = zeros(Nvoxels,oil.ica.num_ics,Nsubs);
group_copes = zeros(Nvoxels,Nics,Ncontrasts);
group_stdcopes = zeros(Nvoxels,Nics,Ncontrasts);
for subnum = 1:Nsubs
    fprintf('Loading in first-level results for subject %d out of %d. \n', subnum, Nsubs);
    copes(:,:,subnum) = nii.quickread(cope_files{subnum},oil.enveloping.gridstep);
end
copes = copes(:,comps2use,:);


wwt_mean = get_weights_norm(oil);
copes = copes./repmat(wwt_mean,[1 size(copes,2),size(copes,3)]);

%% Perform GLM
fprintf('Running GLM. \n');
if demean_designmatrix
    design_matrix=demean(design_matrix,1);
end%if
R       = chol(design_matrix' * design_matrix);
Rinv    = inv(R);
pinvxtx = Rinv * Rinv';
pinvx   = pinvxtx * design_matrix';

for i = 1:Nics
    for v = 1:Nvoxels
        y=permute(copes(v,i,:),[3 1 2]);

        [copeout, varcopeout] = glm_fast_for_meg(y, design_matrix, pinvxtx, ...
                                                 pinvx, contrast,0);
        group_copes(v,i,:)    = copeout;
        group_stdcopes(v,i,:) = sqrt(varcopeout);
    end
end

if group_varcope_spatial_smooth_fwhm > 0
    for c = 1:size(group_stdcopes,3)
        group_stdcopes(:,:,c)=smooth_vol(group_stdcopes(:,:,c), nii.load([OSLDIR '/std_masks/MNI152_T1_' num2str(oil.enveloping.gridstep) 'mm_brain.nii.gz']), group_varcope_spatial_smooth_fwhm, oil.enveloping.gridstep, 'tmp_fname');
    end
end

% Pad out missing components with zeros to allow fslview overlay.
group_tstats = zeros(Nvoxels,oil.ica.num_ics,Ncontrasts);
group_tstats(:,comps2use,:) = group_copes./group_stdcopes;
group_copes_tmp = zeros(Nvoxels,oil.ica.num_ics,Ncontrasts);
group_copes_tmp(:,comps2use,:) = group_copes;
group_copes = group_copes_tmp;

%% Output nifti files
fprintf('Preparing output. \n');
for c = 1:Ncontrasts
    fname_cope = [stats_dir '/copes_group_analysis_components_' num2str(comps2use(1)) '_to_' num2str(comps2use(end)) '_contrast_' num2str(c)];
    fname_tstats=[stats_dir '/tstats_group_analysis_components_' num2str(comps2use(1)) '_to_' num2str(comps2use(end)) '_contrast_' num2str(c)];
    oil.ica_group_level.results.cope_names{c}   = nii.quicksave(group_copes(:,:,c),fname_cope,gridstep,2);
    oil.ica_group_level.results.tstats_names{c} = nii.quicksave(group_tstats(:,:,c),fname_tstats,gridstep,2);
end

%% Permutations testing using RANDOMISE
if use_randomise
    
    fprintf('Permutation testing. \n');
    group_pvals=zeros(Nvoxels,oil.ica.num_ics,Ncontrasts);
    group_clustered_t_stats=zeros(Nvoxels,oil.ica.num_ics,Ncontrasts);
    
    for c=1:Ncontrasts
        dirname = [stats_dir  '/GroupContrast_' num2str(c) '/'];
        if ~isdir(dirname), mkdir(dirname); end
        
        con=contrast{c};
        save_vest(design_matrix,[dirname 'design.mat']);
        save_vest(con',[dirname 'design.con']);
        
        for I = 1:Nics
            fname_cope = [dirname 'copes_component_' num2str(comps2use(I))];
            nii.quicksave(permute(copes(:,I,:),[1 3 2]),fname_cope,oil.enveloping.gridstep);
            output_name = [fname_cope '_output' ];
            
            if ischar(thresh) && strcmp(thresh,'TFC')
            tmp=['randomise -d ' dirname 'design.mat -t ' dirname 'design.con -i ' fname_cope ' -o ' output_name ' -1 -T -R -n ' num2str(nperms) ' --seed=0 -v ' num2str(group_varcope_spatial_smooth_fwhm) ' -m ' mask_fname]; % -c means cluster-based thresholding
            else
            tmp=['randomise -d ' dirname 'design.mat -t ' dirname 'design.con -i ' fname_cope ' -o ' output_name ' -1 -c ' num2str(thresh) ' -R -n ' num2str(nperms) ' --seed=0 -v ' num2str(group_varcope_spatial_smooth_fwhm) ' -m ' mask_fname]; % -c means cluster-based thresholding
            end
            
            disp(tmp);
            runcmd(tmp);
            
            if ischar(thresh) && strcmp(thresh,'TFC')
                group_pvals(:,comps2use(I),c) = nii.quickread([dirname  'copes_component_' num2str(comps2use(I)) '_output_tfce_corrp_tstat1'],gridstep);
                group_clustered_t_stats(:,comps2use(I),c) = nii.quickread([dirname  'copes_component_' num2str(comps2use(I)) '_output_tfce_tstat1'],gridstep);
            else
                group_pvals(:,comps2use(I),c) = nii.quickread([dirname  'copes_component_' num2str(comps2use(I)) '_output_clustere_corrp_tstat1'],gridstep);
                group_clustered_t_stats(:,comps2use(I),c) = nii.quickread([dirname  'copes_component_' num2str(comps2use(I)) '_output_clustere_tstat1'],gridstep);
            end
        end
    end
    
    for c = 1:Ncontrasts
        fname_pvals = [stats_dir '/randomise_corrected_pvals_group_analysis_components_' num2str(comps2use(1)) '_to_' num2str(comps2use(end)) '_contrast_' num2str(c)];
        fname_clusteredtstats=[stats_dir '/randomise_clustered_tstats_group_analysis_components_' num2str(comps2use(1)) '_to_' num2str(comps2use(end)) '_contrast_' num2str(c)];
        oil.ica_group_level.results.corr_1minusp_names{c}   = nii.quicksave(group_pvals(:,:,c),fname_pvals,gridstep,2,'nearestneighbour');
        oil.ica_group_level.results.clustered_tstats_names{c} = nii.quicksave(group_clustered_t_stats(:,:,c),fname_clusteredtstats,gridstep,2,'nearestneighbour');
    end
    
    
    
end

fprintf('Group level OIL complete. \n');
end

function vol_as_matrix=smooth_vol(vol_as_matrix, lower_level_stdbrain, fwhm, gridstep, tmp_fname)

    OSLDIR = getenv('OSLDIR');
    % smooth spatially
    ds=gridstep;

    ss=fwhm/2.3; % std spatial smoothing
    
    nii.save(matrix2vols(vol_as_matrix,lower_level_stdbrain),[ds, ds, ds, 1],[],tmp_fname);

    % Smooth image but respect brain edges.

    maskimage=[OSLDIR '/std_masks/MNI152_T1_' num2str(ds) 'mm_brain_mask'];
    runcmd(['fslmaths ' tmp_fname  ' -s ' num2str(ss) ' -mas ' maskimage ' tmp1']);
    runcmd(['fslmaths ' maskimage ' -s ' num2str(ss)  ' -mas ' maskimage ' tmp2']);
    runcmd(['fslmaths tmp1  -div tmp2 '  tmp_fname]);
    
    runcmd(['rm tmp1.nii* tmp2.nii*']); 

    tmp=nii.load(tmp_fname);
    
    vol_as_matrix=vols2matrix(tmp,lower_level_stdbrain);

end

function wwt_mean = get_weights_norm(oil)
wwt_norm=zeros(3559,numel(oil.concat_subs.sessions_to_do));

for subnum = 1:numel(oil.concat_subs.sessions_to_do)
    
    sing_sub_dat_nn=nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' oil.enveloping.results.source_space_envelopes_NoWeightsNorm_results_fnames{oil.concat_subs.sessions_to_do(subnum)}],oil.enveloping.gridstep);
    sing_sub_dat_wn=nii.quickread([oil.source_recon.dirname '/' oil.enveloping.name '/' oil.enveloping.results.source_space_envelopes_results_fnames{oil.concat_subs.sessions_to_do(subnum)}],oil.enveloping.gridstep);
    
    if subnum==1;
        wwt_norm=zeros(size(sing_sub_dat_nn,1),numel(oil.concat_subs.sessions_to_do));
    end
    wwt_norm(:,subnum)=mean(sing_sub_dat_nn,2)./mean(sing_sub_dat_wn,2);
    
end

wwt_mean = mean(wwt_norm,2);
end
