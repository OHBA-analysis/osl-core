function [ gstats, statsdir ] = osl_parc_permutation_testing( S )
2016-11-22 14:22:30
%%function [ gstats, statsdir ] = osl_parc_permutation_testing( S )
%
% S.oat
% S.time_range % time range
% S.time_average % flag (0 or 1) if 1, this means that the cluster will be
% in 3d, if 0 then will work in 4D.
% S.first_level_copes_to_do=[3]; % list of 1st level contrasts to to perm
% S.group_level_copes_to_do=[1];
% S.parcel_assignments
%

global OSLDIR;

try, masksdir=[OSLDIR '/std_masks' ]; catch, error('OSLDIR not set. Run osl_startup.'); end;

try, S.fsl_version_4p1=S.fsl_version_4p1; catch, S.fsl_version_4p1=1; end;
try, S.matlab_exe_name=S.matlab_exe_name; catch S.matlab_exe_name='matlab'; end;

if(~isfield(S,'time_average'))
    S.time_average=1;
    disp('defaulting to doing time averaging');
end;

statsdir=[];

disp(['Doing cluster perm testing']);

% load in previously run parametric gstats
gstats=oat_load_results(S.oat,S.oat.group_level.results_fnames);

if ~isfield(gstats,'lower_level_copes'),
    warning('Need lower_level_copes to be stored. Re-running group stage to get them');
    oat=S.oat;
    oat.group_level.store_lower_level_copes=1;
    oat.to_do=[0 0 0 1];
    oat=osl_run_oat(oat);
    S.oat=oat;

    gstats=oat_load_results(S.oat,S.oat.group_level.results_fnames);

end;

current_level=S.oat.group_level;

% dummy mask for the fake nifti
nparc = size(gstats.cope,1);
randomise_mask = ones(nparc,1,1);


%% main loop
for coni=1:length(S.first_level_copes_to_do),

    con=S.first_level_copes_to_do(coni);

    cope_smooth_lower_level=gstats.lower_level_copes{con};

    if(size(cope_smooth_lower_level,4)>1)
        error('Not implemented for multiple frequency bins');
    end;

    Sb=[];

    % batch script to run 4D-permutation test on MEG data
    % requires input images for each timepoint stored in a single directory:
    % S.dirname is folder containing a 4D (voxels*subjects) image called
    %       'all_subsXXXX.nii.gz', where XXXX is each timepoint in S.tp

    dirname=[S.oat.source_recon.dirname '/' gstats.fname '_randomise_c' num2str(con) '_dir'];
    mkdir(dirname);

    times=1;
    if isfield(S,'time_range'),
        tinds=intersect(find(gstats.times>S.time_range(1)), find(gstats.times<S.time_range(2)));
        cope_smooth_lower_level=cope_smooth_lower_level(:,:,tinds,:,:);
        if(size(randomise_mask,4)>1),
            randomise_mask=randomise_mask(:,:,:,tinds);
        end;

        times=gstats.times(tinds);
    end;


    if ~S.time_average

        do_tpt=ones(size(cope_smooth_lower_level,3),1);
        for t=1:size(cope_smooth_lower_level,3),

            fnamet=sprintf('%s/allsubs_time%04.0f',dirname,t);
            save_avw(matrix2vols(cope_smooth_lower_level(:,:,t),randomise_mask),fnamet,'f',[nparc, 1, 1, 1]);

            % mask
            fnamet=sprintf('%s/mask_time%04.0f',dirname,t);
            if(size(randomise_mask,4)>1),
                % check mask has any nonzero values at this timepoint
                if ~any(squash(randomise_mask(:,:,:,t)))
                    do_tpt(t)=0;
                else
                    save_avw(randomise_mask(:,:,:,t),fnamet,'f',[nparc, 1, 1, 1]);
                end;
            else
                save_avw(randomise_mask(:,:,:,1),fnamet,'f',[nparc, 1, 1, 1]);
            end;
        end;
    else
        % average over timepoints
        cope_smooth_lower_level=mean(cope_smooth_lower_level,3);
        do_tpt=1;

        fnamet=sprintf('%s/allsubs_time%04.0f',dirname,1);

        save_avw(matrix2vols(cope_smooth_lower_level(:,:,1),randomise_mask),fnamet,'f',[nparc, 1, 1, 1]);

        % mask - use max over timepoints
        disp('Using max of mask over time window');

        randomise_mask=mean(randomise_mask,4);

        fnamet=sprintf('%s/mask_time%04.0f',dirname,1);

        save_avw(randomise_mask(:,:,:,1),fnamet,'f',[nparc, 1, 1, 1]);

    end;

    %  .tp is an array of all timepoints
    Sb.tp=1:size(cope_smooth_lower_level,3);
    % remove timepoints with nothing in the mask
    Sb.tp=Sb.tp(logical(do_tpt));

    %  .nP is number of permutations (must be multiple of 100)
    Sb.nP=S.cluster_stats_nperms;

    %  .thresh is T-statistic threshold to apply to images
    Sb.thresh=S.cluster_stats_thresh;

    Sb.X=current_level.group_design_matrix;

    Sb.gridstep=1;
    Sb.group_varcope_spatial_smooth_std=S.group_varcope_spatial_smooth_fwhm/2.3;

    for gconi=1:length(S.group_level_copes_to_do),
        gcon=S.group_level_copes_to_do(gconi);

        disp(['Cluster 4D perm testing on group contrast ' num2str(gcon)]);
        Sb.contrasts=current_level.group_contrast{gcon};
        Sb.dirname=dirname;

        if S.time_average
            % save design matrix and contrasts
            save_vest(Sb.X',[dirname '/design.mat']);
            save_vest((Sb.contrasts)',[dirname '/design.con']);

            permdir = sprintf('%s',dirname);
            mkdir(permdir);
            tmp=['randomise -d ' dirname '/design.mat -t ' dirname '/design.con -i ' sprintf('%s/allsubs_time%04.0f', Sb.dirname, 1) ' -o ' sprintf('%s/stats', Sb.dirname) ' -R -n ' num2str(Sb.nP) ' --seed=0 -m ' sprintf('%s/mask_time%04.0f', Sb.dirname, 1) ' -x']; % -c means cluster-based thresholding
            disp(tmp);
            runcmd(tmp);

            gstats.dir=Sb.dirname;

            statsdir=permdir;

        else,

            Sb.write_cluster_script=S.write_cluster_script;
            Sb.fsl_version_4p1=S.fsl_version_4p1;
            Sb.times=times;
            Sb.matlab_exe_name=S.matlab_exe_name;
            gstats.clusterstats{con,gcon}=cluster4d_batch(Sb);

            disp('Saving cluster stats.');

            oat_save_results(S.oat,gstats);

            disp('Use osl_save_nii_stats to ouput gstats cluster results.');

        end;

    end;

end

nifs = dir([dirname '/*.nii.gz'])
nii_parcel_settings            = [];
nii_parcel_settings.interp     = 'nearestneighbour';

for idx = 1:length(nifs)
    if strcmp(nifs(idx).name(1:7),'allsubs')
        % Ignore the raw data
        continue
    end
    if strcmp(nifs(idx).name(end-10:end),'parc.nii.gz')
        % Ignore any niftis which are already expanded (and probably about
        % to get overwritten)
        continue
    end

    dat = read_avw([dirname '/' nifs(idx).name]);
    ROInets.nii_parcel_quicksave(dat, S.parcel_assignments, strrep([dirname '/' nifs(idx).name],'.nii.gz','_parc.nii.gz'),nii_parcel_settings);
end


