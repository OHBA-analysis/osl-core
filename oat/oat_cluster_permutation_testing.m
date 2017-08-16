function [ gstats, statsdir ] = oat_cluster_permutation_testing( S )

% [ results ] = oat_cluster_permutation_testing( S )
%
% Runs a group level GLM analysis and uses FSL's randomise to do
% non-parametric stats on the clusters
%
% If S.time_average is 0, then this tests 4D clusters in the timerange specified
% (this is very slow
% and is not run directly, instead it should be only 
% tackled using the scripts that are outputted, on a cluster using
% parallelisation).
%
% If S.time_average is 1, then this tests 3D clusters in the single volume
% in the timerange specified
%
% Input, e.g.:
% S=[];
% S.oat=oat; % oat with a group_level that has been run
% S.time_range % time range
% S.time_average % flag (0 or 1) if 1, this means that the cluster will be
% in 3d, if 0 then will work in 4D.
% S.cluster_stats_thresh=5; % cluster forming threshold
% S.cluster_stats_nperms=5000;
% S.first_level_copes_to_do=[3]; % list of 1st level contrasts to to perm
% S.group_level_copes_to_do=[1];
% testing on
% S.group_varcope_spatial_smooth_fwhm=100; % in mm, this spatially smooths the
% S.write_cluster_script % 0 to run on current processor, or 1 to write a
% script to be run on a cluster (only if 4D)
% group between-subject variances - recommended to be quite high
% S.randomise_mask_fname % mask to limit perms to, can be 4d if perms are
% to be 4d (but can use 3d mask in that context too). 
% Needs to be in same native (low) spatial resolution as oat was run.
%
% MWW 2012

OSLDIR = getenv('OSLDIR');

try, masksdir=[OSLDIR '/std_masks' ]; catch, error('OSLDIR not set. Run osl_startup.'); end;
           
try, S.fsl_version_4p1=S.fsl_version_4p1; catch, S.fsl_version_4p1=1; end;
try, S.matlab_exe_name=S.matlab_exe_name; catch S.matlab_exe_name='matlab'; end;

if(~isfield(S,'time_average'))
    S.time_average=1;
    disp('defaulting to doing time averaging');
end;

if(isfield(S,'timepoint'))
    error('S.timepoint is a deprecated option. Use S.time_range and S.time_average instead');
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

% sort out oat group masks (in native res)
% current_level_mask_fname=[S.oat.source_recon.dirname '/' current_level.name '_mask'];
current_level_mask_fname=[S.oat.source_recon.dirname '/' S.oat.first_level.name '_' S.oat.subject_level.name '_' current_level.name '_mask.nii.gz'];
[std_brainmask,std_res,std_xform]=nii.load(current_level_mask_fname);

% load in mask to use for randomise (assumed to be in native res)
if isfield(S,'randomise_mask_fname'),
    [randomise_mask,randomise_res,randomise_xform]=nii.load(S.randomise_mask_fname);
else
    randomise_mask=std_brainmask;
    randomise_res = std_res;
    randomise_xform = std_xform;
end;

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
            nii.save(matrix2vols(cope_smooth_lower_level(:,:,t),std_brainmask),std_res,std_xform,fnamet);

            % mask
            fnamet=sprintf('%s/mask_time%04.0f',dirname,t);
            if(size(randomise_mask,4)>1),        
                % check mask has any nonzero values at this timepoint
                if ~any(squash(randomise_mask(:,:,:,t)))
                    do_tpt(t)=0;  
                else
                    nii.save(randomise_mask(:,:,:,t),randomise_res,randomise_xform,fnamet);
                end;
            else
                nii.save(randomise_mask(:,:,:,1),randomise_res,randomise_xform,fnamet);
            end;
        end;
    else
        % average over timepoints
        cope_smooth_lower_level=mean(cope_smooth_lower_level,3);
        do_tpt=1;
        
        fnamet=sprintf('%s/allsubs_time%04.0f',dirname,1);
        
        nii.save(matrix2vols(cope_smooth_lower_level(:,:,1),std_brainmask),std_res,std_xform,fnamet);

        % mask - use max over timepoints
        disp('Using max of mask over time window');
        
        randomise_mask=mean(randomise_mask,4);
        
        fnamet=sprintf('%s/mask_time%04.0f',dirname,1);
        
        nii.save(randomise_mask(:,:,:,1),randomise_res,randomise_xform,fnamet);
         
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

    Sb.gridstep=gstats.gridstep;
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
            tmp=['randomise -d ' dirname '/design.mat -t ' dirname '/design.con -i ' sprintf('%s/allsubs_time%04.0f', Sb.dirname, 1) ' -o ' sprintf('%s/stats', Sb.dirname) ' -c ' num2str(Sb.thresh) ' -R -n ' num2str(Sb.nP) ' --seed=0 -v ' num2str(Sb.group_varcope_spatial_smooth_std) ' -m ' sprintf('%s/mask_time%04.0f', Sb.dirname, 1)]; % -c means cluster-based thresholding
            disp(tmp);
            runcmd(tmp);
            
            gstats.dir=Sb.dirname;
            
            gridstep=gstats.gridstep;
            
            resamp_gridstep=gridstep;
  
            origname='tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','cubic','enforce_mask',true);
            origname='clustere_tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','nearest','enforce_mask',true);
            origname='clustere_corrp_tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','nearest','enforce_mask',true);
      
            resamp_gridstep=2;
            
            origname='tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','cubic','enforce_mask',true);
            origname='clustere_tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','nearest','enforce_mask',true);
            origname='clustere_corrp_tstat';
            nii.resample([permdir '/stats_' origname num2str(1) '.nii.gz'],[permdir '/' origname num2str(con) '_gc' num2str(gcon) '_' num2str(resamp_gridstep) 'mm.nii.gz'],resamp_gridstep,'interptype','nearest','enforce_mask',true);
            
            statsdir=permdir;
            
        else,
            
            Sb.write_cluster_script=S.write_cluster_script;
            Sb.fsl_version_4p1=S.fsl_version_4p1;
            Sb.times=times;
            Sb.matlab_exe_name=S.matlab_exe_name;
            gstats.clusterstats{con,gcon}=cluster4d_batch(Sb);
            
            disp('Saving cluster stats.');
            
            oat_save_results(S.oat,gstats);

            disp('Use oat_save_nii_stats to ouput gstats cluster results.');
            
        end;
    end;

end;



