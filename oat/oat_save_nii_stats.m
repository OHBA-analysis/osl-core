function [statsdir,times,cluster_times]=oat_save_nii_stats( Sin )

% [statsdir,times,cluster_times]=oat_save_nii_stats( Sin )
%
% Outputs nii files for the specified first or group level contrasts from an OAT analysis 
% into statsdir. Can be used for results from other ( none OAT ) analyses as well.
%
% By default nii files are output in the native and 2mm resolution. Can be used for first or group level
% stats. For group level all group level contrasts are ouput for the
% specified first level contrasts.
%
% Required inputs:
% Sin.oat: oat struct resulting from an OAT analysis
%
% Sin.stats_fname: first_level or group_level results struct file name, e.g. resulting from
%                   an OAT analysis
% OR
% Sin.stats: first_level or group_level results struct, e.g. resulting from
%               an OAT analysis
%
% Sin.stats_dir: specify ouput directory to put nii files into
% Sin.stdbrainmask_fname: std brain mask at grid resolution of the analysis
%
% Optional inputs:
% Sin.stats_dir: specify ouput directory to put nii files into
% Sin.first_level_contrasts: list of first level contrast indices
%                               to indicate which first level contrasts to output
% Sin.group_level_contrasts: list of group level contrast indices
%                               to indicate which group level contrasts to output
%                               should have nii files output for them
% Sin.resamp_gridstep: resolution to resample niis to [default is 2mm]
% Sin.time_range: time range to use, i.e. [from to]. Default is to use full
%                   time range.
% Sin.freq_bin: index of freq bin to us [default is 1]
%
% Example first level stats call:
%
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.first_level.results_fnames{1};
% S2.first_level_contrasts=[1,2]; % list of first level contrasts to output
% statsdir=oat_save_nii_stats(S2);
%
% Example group level stats call:
% 
% S2=[];
% S2.oat=oat;
% S2.stats_fname=oat.group_level.results_fnames;
% S2.first_level_contrasts=[1,2]; % list of first level contrasts to output
% statsdir=oat_save_nii_stats(S2);
%
% Example none OAT call:
%
% S2=[];
% S2.stdbrainmask_fname=stdbrainmask_fname;
% S2.stats_dir=stats_dir;
% stats.cope=cope;
% stats.times=times;
% S2.stats=stats; % list of first level contrasts to output
% statsdir=oat_save_nii_stats(S2);
%
% where cope is num_voxels x num_timepoints
%
% MWW 2011

global OSLDIR;

try, resamp_gridstep=Sin.resamp_gridstep; catch, resamp_gridstep=2; end; % mm

try, masksdir=[OSLDIR '/std_masks' ]; catch, error('OSLDIR not set. Run osl_startup.'); end;

try, resolve_erf_sign_ambiguity=Sin.resolve_erf_sign_ambiguity; catch, Sin.resolve_erf_sign_ambiguity=0; end;

try, tmp=Sin.freq_bin; catch, Sin.freq_bin=1; end;

try, Sin.output_varcopes=Sin.output_varcopes; catch, Sin.output_varcopes=0; end;

Sin.have_stdcope=1;

cluster_times=[];

if(isfield(Sin,'oat')),
    if(isfield(Sin,'stats_fname')),
        Sin.stats=oat_load_results(Sin.oat,Sin.stats_fname);            
    end;

    try, Sin.stats_dir=Sin.stats_dir; catch, 
        Sin.stats_dir=[Sin.oat.source_recon.dirname '/' Sin.stats.fname '_dir']; 
    end; %full path of directory where stats will be stored
   
    if Sin.stats.level==1,
        stdbrainmask_fname=[Sin.oat.source_recon.dirname '/' Sin.oat.first_level.name '_mask'];
    elseif Sin.stats.level==2,
        stdbrainmask_fname=[Sin.oat.source_recon.dirname '/' Sin.oat.first_level.name '_' Sin.oat.subject_level.name '_' Sin.oat.group_level.name '_mask'];
    else
        error('Invalid stats.level');
    end;
else

    try, Sin.stats_dir=Sin.stats_dir; catch, error('Sin.stats_dir not specified'); end; %full path of directory where stats will be stored
    try, stdbrainmask_fname=Sin.stdbrainmask_fname; catch, error('Sin.stdbrainmask_fname not specified'); end;
    
    if(~isfield(Sin.stats,'level')),
        Sin.stats.level=1;
    end;
    
    if(~isfield(Sin.stats,'stdcope')),
        Sin.have_stdcope=0;
    end;
end;

try, Sin.first_level_contrasts=Sin.first_level_contrasts; catch, first_level_contrasts=1:size(Sin.stats.cope,3); end; 
try, Sin.group_level_contrasts=Sin.group_level_contrasts; catch, group_level_contrasts=1:size(Sin.stats.cope,5); end; 

if(isfield(Sin.stats,'clusterstats')),
    try,corrp_thresh=Sin.cluster_stats_fpr;catch, corrp_thresh=0.05;end;
end;

tim=Sin.stats.times;

if(length(tim)>1)
    tres=tim(2)-tim(1);
else
    tres=1;
end;

if(Sin.resolve_erf_sign_ambiguity),
    Sin.stats=osl_resolve_sign_ambiguity(Sin);
    count=Sin.stats.count;
else
    count=0;
end;

if(isfield(Sin,'time_range')),
    time_indices=find(tim>=Sin.time_range(1) & tim<=Sin.time_range(2));
else
    time_indices=1:length(tim);   
end;

clear tim;

if(size(Sin.stats.cope,4)>1 && ~isfield(Sin,'freq_bin'))
    error('Not implemented for multiple freq bins');
end;

warning off;
mkdir(Sin.stats_dir);
warning on;

% setup options for nii_quicksave calls
options=[];
options.interp='trilinear';
options.mask_fname=stdbrainmask_fname;
options.output_spat_res=resamp_gridstep;
options.tres=tres;

nii_parcel_settings            = [];
nii_parcel_settings.interp     = 'nearestneighbour';
%nii_parcel_settings.mask_fname = Sin.stats.D_sensor_data.parcellation.S.parcellation;
use_parcels=0;

try
    use_parcels=isfield(Sin.stats.D_sensor_data,'parcellation');
catch
end
      
if(size(Sin.stats.cope,4)>1)
    Sin.stats.cope=Sin.stats.cope(:,:,:,Sin.freq_bin,:);
    Sin.stats.stdcope=Sin.stats.stdcope(:,:,:,Sin.freq_bin,:);
end;

if(Sin.stats.level==1),

    if(isfield(Sin.stats,'pseudo_zstat_var')),
        fname=[Sin.stats_dir '/pseudo_zstat_var'];
        dat=Sin.stats.pseudo_zstat_var;
        
        if use_parcels                
            fname_out = ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
        else
            fname_out = nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);
        end
    end;
    
    if(isfield(Sin.stats,'sqrt_wnorm_nai')),
        fname=[Sin.stats_dir '/sqrt_wnorm_nai'];    
        dat=Sin.stats.sqrt_wnorm_nai;
        if use_parcels                
            fname_out = ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
        else
            fname_out = nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);
        end
    end;
    
    gcon=1;
    fname_postfix='';    
    fnamet=output_files(Sin,options,time_indices,gcon,fname_postfix,use_parcels,nii_parcel_settings);
    
elseif(Sin.stats.level==2)     
    
    for gconi=1:length(group_level_contrasts),        
        gcon=group_level_contrasts(gconi);
        fname_postfix=['_gc' num2str(gcon)];        
        fnamet=output_files(Sin,options,time_indices,gcon,fname_postfix,use_parcels,nii_parcel_settings);        
    end;  
    
else
    
    error('level unsupported');
    
end

times=Sin.stats.times(time_indices);
save([Sin.stats_dir '/times'],'times','-ascii');
statsdir=Sin.stats_dir;

disp(['E.g. to view nii file: ']);
disp(['fslview(''' fnamet ''')']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnamet=output_files(Sin,options,time_indices,gcon,fname_postfix,use_parcels,nii_parcel_settings)

statsdir=Sin.stats_dir;
for coni=1:length(Sin.first_level_contrasts),

    con=Sin.first_level_contrasts(coni);

    disp(['Outputting nii files for first level contrast ' num2str(con)]);    

    % output nifti file for the contrast of parameter estimates
    % (COPEs) for each contrast
    fname=[statsdir '/cope' num2str(con) fname_postfix];
    dat=permute(Sin.stats.cope(:,time_indices,con,:,gcon),[1 2 4 3]);
    
    if use_parcels                
        ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
    else
        nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);        
    end
    
    fnamet=fname;
    
    if(Sin.have_stdcope)      

        fname=[statsdir '/tstat' num2str(con) fname_postfix];
        ts=Sin.stats.cope(:,time_indices,con,:,gcon)./Sin.stats.stdcope(:,time_indices,con,:,gcon);
        dat=permute(ts,[1 2 4 3]);
        
        if use_parcels                
            fnamet=ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
        else
            fnamet=nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);        
        end
        
        if Sin.output_varcopes,
            fname=[statsdir '/varcope' num2str(con) fname_postfix]; 
            dat=permute(Sin.stats.stdcope(:,time_indices,con,:,gcon).^2,[1 2 4 3]);
            if use_parcels                
                fnamet=ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
            else
                fnamet=nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);
            end
        end;

        if size(ts,2)>1,

            fname=[statsdir '/tstat' num2str(con) '_mip' fname_postfix];
            times=Sin.stats.times(time_indices);
            time_indices_sum=find(times>=0);
            time_indices_sum=intersect(time_indices_sum,time_indices);
            ts=abs(Sin.stats.cope(:,time_indices_sum,con,:,gcon)./Sin.stats.stdcope(:,time_indices_sum,con,:,gcon));

            % output max
            ts=max(ts,[],2);

            dat=permute(ts,[1 2 4 3]);
            
            if use_parcels                
                fnamet_mip=ROInets.nii_parcel_quicksave(dat, Sin.stats.D_sensor_data.parcellation.assignments, [fname '_' num2str(options.output_spat_res) 'mm'], nii_parcel_settings);
            else
                fnamet_mip=nii_quicksave(dat,[fname '_' num2str(options.output_spat_res) 'mm'],options);
            end

        end;
    end;
    
    if(isfield(Sin.stats,'clusterstats')),
                
        disp('Outputting nii files for cluster stats');

        nC = length(Sin.stats.clusterstats{con,gcon}.Creal);
        pVreal = zeros(nC,1);
        pVimg = Sin.stats.clusterstats{con,gcon}.clustimg;
        tstat_thresh_img = zeros(size(Sin.stats.clusterstats{con,gcon}.clustimg));

        for i = 1:nC %loop over real clusters

            % Sin.stats.clusterstats{con,gcon}.nVreal is a list of the #voxels in each cluster 
            % Sin.stats.clusterstats{con,gcon}.dist is the null distribution of cluster extent
            % Sin.stats.clusterstats{con,gcon}.Creal is a list of cluster index for each cluster
            % Sin.stats.clusterstats{con,gcon}.clustimg is a 4D image of cluster indices
            pVreal(i) = mean(Sin.stats.clusterstats{con,gcon}.nVreal(i)>=Sin.stats.clusterstats{con,gcon}.dist); % this is 1-P-Value for each cluster
            pVimg(Sin.stats.clusterstats{con,gcon}.clustimg==full(Sin.stats.clusterstats{con,gcon}.Creal(i))) = full(pVreal(i));

        end;

        clustimg_fname = [fnamet '_clust4d'];
        pVimg_fname = [fnamet '_clust4d_corrp'];
        tstats_clust_fname = [fnamet '_clust4d_tstats'];
        cluster_times=Sin.stats.clusterstats{con,gcon}.times;

        save_avw(Sin.stats.clusterstats{con,gcon}.clustimg,[clustimg_fname '_' num2str(gridstep) 'mm'],'i',[gridstep gridstep gridstep tres]);
        save_avw(pVimg,[pVimg_fname '_' num2str(gridstep) 'mm'],'f',[gridstep gridstep gridstep tres]);
        save_avw(Sin.stats.clusterstats{con,gcon}.tstats,[tstats_clust_fname '_' num2str(gridstep) 'mm'],'f',[gridstep gridstep gridstep tres]);

        osl_resample_nii([clustimg_fname '_' num2str(gridstep) 'mm'],[clustimg_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'trilinear',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);                    
        osl_resample_nii([pVimg_fname '_' num2str(gridstep) 'mm'],[pVimg_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'nearestneighbour',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);
        osl_resample_nii([tstats_clust_fname '_' num2str(gridstep) 'mm'],[tstats_clust_fname '_' num2str(resamp_gridstep) 'mm'],resamp_gridstep,'nearestneighbour',[masksdir '/MNI152_T1_' num2str(resamp_gridstep) 'mm_brain_mask' ],0);
        runcmd(['fslmaths ' clustimg_fname '_' num2str(resamp_gridstep) 'mm -thr ' num2str(1-corrp_thresh) ' ' clustimg_fname '_thresh_' num2str(resamp_gridstep) 'mm']);
        % runcmd(['fslmaths ' fnamet '_' num2str(resamp_gridstep) 'mm -mas ' clustimg_fname '_thresh_' num2str(resamp_gridstep) 'mm -thr ' num2str(Sin.stats.S.cluster_stats_thresh) ' ' tstat_thresh_img_fname '_' num2str(resamp_gridstep) 'mm ']);
    end
   
end